import os
import sys
from sys import argv
import argparse
import logging
import subprocess
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger


def concatenate_hits(logger, normalized_hits_dir, job_input_file, output_dir):
    job_index = os.path.splitext(os.path.basename(job_input_file))[0]
    output_file_path = os.path.join(output_dir, f"temp_{job_index}.txt")
    # avoid cat {input_dir}/* because arguments list might be too long!
    # No need to wait...

    with open(job_input_file, 'r') as f:
        for line in f:
            hits_file_name = line.strip()
            hits_file_path = os.path.join(normalized_hits_dir, hits_file_name)

            cmd = f"cat {hits_file_path} >> {output_file_path}"
            logger.info(f'Calling: {cmd}')
            subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('normalized_hits_dir', help='dir of hits to concatenate')
    parser.add_argument('job_input_file', help='path to a file that contains the names of the files to concatenate')
    parser.add_argument('output_dir', help='output dir of concatenated files')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        concatenate_hits(logger, args.normalized_hits_dir, args.job_input_file, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
