import pandas as pd
from sys import argv
import argparse
import os
import sys
import logging
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def normalize_hits_scores(logger, blast_result, output_path, scores_normalize_coefficient, use_parquet):
    logger.info(f'Normalize scores of {blast_result}')

    query_vs_reference_file_name = os.path.splitext(os.path.basename(blast_result))[0]
    strain1_name, strain2_name = query_vs_reference_file_name.split('_vs_')

    if use_parquet:
        df = pd.read_parquet(blast_result)
    else:
        df = pd.read_csv(blast_result)

    df['score'] = df['score'] / scores_normalize_coefficient
    df.to_csv(output_path, index=False, header=[strain1_name, strain2_name, 'score'])


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_result', help='path to filtered mmseqs result')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('scores_normalize_coefficient')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--use_parquet', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        normalize_hits_scores(logger, args.blast_result, args.output_path,
                              float(args.scores_normalize_coefficient), args.use_parquet)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
