import argparse
import sys
import time
from math import floor
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.general_utils import remove_path
from pipeline.flask import SharedConsts, flask_interface_consts


def clean_old_jobs(logger):
    if not flask_interface_consts.USER_RESULTS_DIR.exists():
        logger.info(f'{flask_interface_consts.USER_RESULTS_DIR} is not an existing directory.')
        return

    for job_dir_path in flask_interface_consts.USER_RESULTS_DIR.iterdir():
        if not job_dir_path.is_dir():
            continue

        seconds_since_last_modification = time.time() - job_dir_path.stat().st_mtime
        days_since_last_modification = floor(seconds_since_last_modification / 60 / 60 / 24)

        if days_since_last_modification > SharedConsts.TIME_TO_KEEP_PROCSES_IDS_FOLDERS + 1:  # Take confidence interval of 1 day
            logger.info(
                f'Removing old job dir: {job_dir_path}, since it was done {days_since_last_modification} days ago')
            remove_path(logger, job_dir_path)
            logger.info(f'Removed {job_dir_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, clean_old_jobs)
