"""
Weekly job report for Microbializer production server.

Scans USER_RESULTS_DIR for jobs from the last 7 days, classifies each job as
finished / bad_input / internal_error, and emails a summary to the owners.

Error classification (based on error.txt content):
  - No error.txt + finished marker exists  -> finished
  - error.txt contains a Python traceback
    or "Error in command:" prefix           -> internal_error (needs investigation)
  - error.txt exists with only readable text -> bad_input (user's fault)
  - Neither file exists                      -> still running (skipped)
"""

import re
import sys
import logging
import argparse
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Make pipeline importable when run directly
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.flask.flask_interface_consts import (
    USER_RESULTS_DIR,
    ERROR_FILE_PATH,
    FINISHED_JOB_FILE_PATH,
    JOB_PARAMETERS_FILE_NAME,
    ADDITIONAL_OWNER_EMAILS,
)
from pipeline.flask.SharedConsts import WEBSERVER_ADDRESS, OWNER_EMAIL
from pipeline.auxiliaries.email_sender import send_email

OWNER_EMAILS = [OWNER_EMAIL] + ADDITIONAL_OWNER_EMAILS

# Patterns that indicate a non-user-facing (internal) error
_INTERNAL_ERROR_PATTERNS = re.compile(
    r'Traceback \(most recent call last\):|Error in command:',
    re.MULTILINE,
)


def _job_start_time(job_dir: Path) -> datetime:
    """Best estimate of when the job was submitted."""
    params_file = job_dir / JOB_PARAMETERS_FILE_NAME
    ref = params_file if params_file.exists() else job_dir
    return datetime.fromtimestamp(ref.stat().st_mtime)


def _classify_job(job_dir: Path) -> Optional[str]:
    """
    Returns 'finished', 'bad_input', 'internal_error', or None (still running).
    """
    finished_marker = job_dir / FINISHED_JOB_FILE_PATH
    error_file = job_dir / ERROR_FILE_PATH

    if finished_marker.exists():
        return 'finished'

    if error_file.exists():
        content = error_file.read_text(errors='replace')
        if _INTERNAL_ERROR_PATTERNS.search(content):
            return 'internal_error'
        return 'bad_input'

    return None


def _build_email_body(
    jobs_by_status: Dict[str, List[Tuple[str, datetime]]],
    week_start: datetime,
    now: datetime,
) -> str:
    finished = jobs_by_status['finished']
    bad_input = jobs_by_status['bad_input']
    internal_error = jobs_by_status['internal_error']
    total = len(finished) + len(bad_input) + len(internal_error)

    lines = [
        'Microbializer Weekly Job Report',
        f'Period: {week_start.strftime("%Y-%m-%d %H:%M")} — {now.strftime("%Y-%m-%d %H:%M")}',
        '',
        'Summary',
        '-------',
        f'  Total jobs submitted:  {total}',
        f'  Finished successfully: {len(finished)}',
        f'  Bad input (user error): {len(bad_input)}',
        f'  Internal error (action required): {len(internal_error)}',
    ]

    def _section(title, jobs, url_path):
        if not jobs:
            return
        lines.append('')
        lines.append(f'=== {title} ({len(jobs)}) ===')
        for job_id, start_time in jobs:
            lines.append(
                f'  {start_time.strftime("%Y-%m-%d %H:%M")}  '
                f'{WEBSERVER_ADDRESS}/{url_path}/{job_id}'
            )

    _section('Finished Jobs', finished, 'results')
    _section('Bad Input Errors — no action needed', bad_input, 'error_from_job')
    _section('Internal Errors — REQUIRES INVESTIGATION', internal_error, 'error_from_job')

    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Send Microbializer weekly job report.')
    parser.add_argument(
        '--days', type=int, default=7,
        help='Number of past days to include in the report (default: 7).',
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Print the email body to stdout instead of sending it.',
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    logger = logging.getLogger('weekly_report')

    now = datetime.now()
    cutoff = now - timedelta(days=args.days)

    user_results_dir = Path(USER_RESULTS_DIR)
    if not user_results_dir.exists():
        logger.error('USER_RESULTS_DIR does not exist: %s', user_results_dir)
        sys.exit(1)

    jobs_by_status: Dict[str, List[Tuple[str, datetime]]] = {
        'finished': [],
        'bad_input': [],
        'internal_error': [],
    }

    for job_dir in user_results_dir.iterdir():
        if not job_dir.is_dir():
            continue

        start_time = _job_start_time(job_dir)
        if start_time < cutoff:
            continue

        status = _classify_job(job_dir)
        if status is None:
            logger.info('Skipping still-running job: %s', job_dir.name)
            continue

        jobs_by_status[status].append((job_dir.name, start_time))

    for status in jobs_by_status:
        jobs_by_status[status].sort(key=lambda x: x[1])

    subject = (
        f'Microbializer Weekly Report — {now.strftime("%Y-%m-%d")} '
        f'({len(jobs_by_status["internal_error"])} internal errors)'
    )
    body = _build_email_body(jobs_by_status, cutoff, now)

    if args.dry_run:
        print(f'Subject: {subject}')
        print()
        print(body)
        return

    logger.info('Sending report to %s', OWNER_EMAILS)
    send_email(logger, subject, body, OWNER_EMAILS)
    logger.info('Done.')


if __name__ == '__main__':
    main()
