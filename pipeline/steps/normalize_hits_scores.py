import pandas as pd
from sys import argv
import argparse
import os
import sys
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts


def normalize_hits_scores(logger, filtered_blast_result, output_path, scores_normalize_coefficient):
    logger.info(f'Normalize scores of {filtered_blast_result}')

    df = pd.read_csv(filtered_blast_result)
    df['score'] = df['score'] / scores_normalize_coefficient
    df.to_csv(output_path, index=False)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('filtered_blast_result', help='path to filtered mmseqs result')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('scores_normalize_coefficient')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        normalize_hits_scores(logger, args.filtered_blast_result, args.output_path,
                              float(args.scores_normalize_coefficient))
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
