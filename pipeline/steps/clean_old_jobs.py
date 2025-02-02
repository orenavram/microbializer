from sys import argv
import argparse
import sys
import traceback
import os
import time
from math import floor
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries import consts
from auxiliaries.pipeline_auxiliaries import get_job_logger, remove_path, add_default_step_args
from flask import SharedConsts


def clean_old_jobs(logger):
    if not consts.USER_RESULTS_DIR.exists():
        logger.info(f'{consts.USER_RESULTS_DIR} is not an existing directory.')
        return

    for job_dir_path in consts.USER_RESULTS_DIR.iterdir():
        if not job_dir_path.isdir():
            continue

        seconds_since_last_modification = time.time() - os.stat(job_dir_path).st_mtime
        days_since_last_modification = floor(seconds_since_last_modification / 60 / 60 / 24)

        if days_since_last_modification > SharedConsts.TIME_TO_KEEP_PROCSES_IDS_FOLDERS + 1:  # Take confidence interval of 1 day
            logger.info(f'Removing old job dir: {job_dir_path}, since it was done {days_since_last_modification} days ago')
            remove_path(logger, job_dir_path)
            logger.info(f'Removed {job_dir_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        clean_old_jobs(logger)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
