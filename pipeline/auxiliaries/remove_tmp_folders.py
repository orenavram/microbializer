from sys import argv
import argparse
import shutil
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def remove_dirs_from_tmp_dir(logger, tmp_dir):
    for file in os.listdir(tmp_dir):
        folder = os.path.join(tmp_dir, file)
        if os.path.isdir(folder):
            logger.info(f'Removing {folder}')
            shutil.rmtree(folder, ignore_errors=True)


if __name__ == '__main__':
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    parser = argparse.ArgumentParser()
    parser.add_argument('tmp_dir', help='A path from which all folders will be removed')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir)

    remove_dirs_from_tmp_dir(logger, args.tmp_dir)
