from sys import argv
import argparse
import logging
import os
import sys
import traceback
import time
from math import floor

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries import consts
from auxiliaries.pipeline_auxiliaries import get_job_logger, remove_path
from flask import SharedConsts


def clean_old_jobs(logger):
    if not os.path.isdir(consts.USER_RESULTS_DIR):
        logger.info(f'{consts.USER_RESULTS_DIR} is not an existing directory.')
        return

    for job_dir_name in os.listdir(consts.USER_RESULTS_DIR):
        job_dir_path = os.path.join(consts.USER_RESULTS_DIR, job_dir_name)
        if not os.path.isdir(job_dir_path):
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
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        clean_old_jobs(logger)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
