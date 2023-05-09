from sys import argv
import argparse
import logging
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def verify(logger, input_file, output_file, clustering_criterion):
    with open(input_file) as f:
        i = 0
        for line in f:
            # each line is a single cluster
            i += 1
            if i > clustering_criterion:
                logger.info(f'{input_file} has at least {i} clusters, thus not relevant.')
                return False
        if i == 0:
            raise ValueError(f'{input_file} is empty! There\'s a bug in the previous step!')

    os.rename(input_file, output_file)
    return True


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('input_file', help='path to an MCL analysis file')
    parser.add_argument('output_file',
                        help='path to which the MCL analysis will be moved if clustering criterion was met')
    parser.add_argument('--clustering-criterion', help='maximal number of clusters allowed', type=int, default=1)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        verify(logger, args.input_file, args.output_file, args.clustering_criterion)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
