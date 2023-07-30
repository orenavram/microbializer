from sys import argv
import argparse
import logging
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def verify(logger, input_file, output_dir, clustering_criterion):
    input_og_name = os.path.splitext(os.path.basename(input_file))[0]
    with open(input_file, 'r') as f:
        lines = [line.rstrip() for line in f]

    if len(lines) > clustering_criterion:
        logger.info(f'{input_file} has {len(lines)} clusters, thus not relevant.')
    elif len(lines) == 0:
        raise ValueError(f'{input_file} is empty! There\'s a bug in the previous step!')
    else:  # 1 <= len(lines) <= clustering_criterion
        for line in lines:
            genes = line.split('\t')
            if input_og_name in genes:
                og_name = input_og_name
            else:
                og_name = genes[0]
            verified_cluster_path = os.path.join(output_dir, og_name + ".10_verified_cluster")
            with open(verified_cluster_path, 'w') as verified_cluster_file:
                verified_cluster_file.write(line)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='path to an MCL analysis file')
    parser.add_argument('output_dir',
                        help='dir path to which the MCL analysis will be moved if clustering criterion was met')
    parser.add_argument('--clustering-criterion', help='maximal number of clusters allowed', type=int, default=4)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        verify(logger, args.input_file, args.output_dir, args.clustering_criterion)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
