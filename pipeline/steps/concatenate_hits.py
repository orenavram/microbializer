import os
import sys
from sys import argv
import argparse
import logging
import subprocess
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
from auxiliaries.logic_auxiliaries import max_with_nan


def concatenate_hits(logger, input_dir, start_index, end_index, output_dir):
    output_file_path = os.path.join(output_dir, f"temp_{start_index}_{end_index}.txt")
    # avoid cat {input_dir}/* because arguments list might be too long!
    # No need to wait...
    for file in os.listdir(input_dir)[start_index:end_index]:
        cmd = f"cat {input_dir}/{file} >> {output_file_path}"
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir', help='dir of hits to concatenate')
    parser.add_argument('start_index', help='start index of files list', type=int)
    parser.add_argument('end_index', help='end index of files list', type=int)
    parser.add_argument('output_dir', help='output dir of concatenated files')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        concatenate_hits(logger, args.input_dir, args.start_index, args.end_index, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')