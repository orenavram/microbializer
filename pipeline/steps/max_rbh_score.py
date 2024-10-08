import os
import sys
from sys import argv
import argparse
import logging
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
from auxiliaries.logic_auxiliaries import max_with_nan


def max_rbh_score_per_gene(logger, rbh_m8_dir, strain_name, output_dir, step_name, error_file_path):
    max_score_per_gene = pd.Series(dtype=float)

    output_file_path = os.path.join(output_dir, f'{strain_name}.{step_name}')
    if not os.path.exists(output_file_path):
        for rbh_hits_file in os.listdir(rbh_m8_dir):
            if 'm8' not in rbh_hits_file:
                continue
            query_vs_reference_file_name = os.path.splitext(rbh_hits_file)[0]
            query_strain, target_strain = query_vs_reference_file_name.split('_vs_')
            if query_strain != strain_name and target_strain != strain_name:
                continue
            try:
                rbh_hits_df = pd.read_csv(os.path.join(rbh_m8_dir, rbh_hits_file), sep='\t')
                if query_strain == strain_name:
                    genes_max_scores = rbh_hits_df.groupby(['query']).max(numeric_only=True)['score']
                else: # target_strain == strain_name
                    genes_max_scores = rbh_hits_df.groupby(['target']).max(numeric_only=True)['score']

                max_score_per_gene = max_score_per_gene.combine(genes_max_scores, max_with_nan)
            except Exception as e:
                logger.exception(f'Error while processing {rbh_hits_file} in step max_rbh_score_per_gene: {e}')
                fail(logger, f'Error while processing {rbh_hits_file} in step max_rbh_score_per_gene: {e}', error_file_path)

        max_score_per_gene.to_csv(output_file_path, index_label='gene', header=['max_ortholog_score'])


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('rbh_m8_dir', help='dir of m8 files of rbh analysis')
    parser.add_argument('strain_name', help='strain name')
    parser.add_argument('output_dir', help='path to which the results will be written')
    parser.add_argument('step_name', help='step name')
    parser.add_argument('error_file_path', help='path to which errors are written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        max_rbh_score_per_gene(logger, args.rbh_m8_dir, args.strain_name, args.output_dir, args.step_name, args.error_file_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            f.write(f'Internal Error in {__file__}: {e}\n')
