from sys import argv
import argparse
import logging
import subprocess
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def mcl(logger, input_file, output_file):
    # --abc for a columns format, i.e., item1\item2\tscore
    cmd = f'mcl "{input_file}" --abc -o "{output_file}"'
    logger.info(f'Starting MCL. Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('input_file', help='path to an MCL input file')
    parser.add_argument('output_file', help='path to which the MCL analysis will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    mcl(logger, args.input_file, args.output_file)
