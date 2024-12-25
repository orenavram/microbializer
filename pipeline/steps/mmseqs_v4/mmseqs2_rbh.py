import os
import subprocess
import sys
import time
from sys import argv
import argparse
import logging
import shutil
import pandas as pd
import traceback
import statistics
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
from auxiliaries.logic_auxiliaries import add_score_column_to_mmseqs_output
from auxiliaries import consts


def search_rbh(logger, genome1, genome2, dbs_dir, rbh_hits_dir, scores_statistics_dir, max_rbh_score_per_gene_dir,
               temp_dir, identity_cutoff, coverage_cutoff, e_value_cutoff, use_parquet):
    """
    input:  2 protein dbs
    output: query_vs_reference "mmseqs2 easy-rbh" results file
    """
    output_rbh_path = os.path.join(rbh_hits_dir, f'{genome1}_vs_{genome2}.m8')
    output_statistics_path = os.path.join(scores_statistics_dir, f'{genome1}_vs_{genome2}.stats')
    output_genome1_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome1}_max_scores_with_{genome2}.csv')
    output_genome2_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome2}_max_scores_with_{genome1}.csv')

    if os.path.exists(output_rbh_path) and os.path.exists(output_statistics_path) \
            and os.path.exists(output_genome1_max_scores) and os.path.exists(output_genome2_max_scores):
        return

    tmp_dir = os.path.join(temp_dir, f'tmp_{genome1}_vs_{genome2}')

    # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
    db_1_path = os.path.join(dbs_dir, f'{genome1}.db')
    db_2_path = os.path.join(dbs_dir, f'{genome2}.db')
    result_db_path = os.path.join(temp_dir, f'{genome1}_vs_{genome2}.db')
    rbh_command = f'mmseqs rbh {db_1_path} {db_2_path} {result_db_path} {tmp_dir} --min-seq-id {identity_cutoff} ' \
                  f'-c {coverage_cutoff} --cov-mode 0 -e {e_value_cutoff} --threads 1 --search-type 1 ' \
                  f'--comp-bias-corr 0 -v 1 -s {consts.MMSEQS_SENSITIVITY_PARAMETER}'
    logger.info(f'Calling: {rbh_command}')
    subprocess.run(rbh_command, shell=True, check=True)

    if os.path.getsize(result_db_path) == 0:
        logger.info(f"{result_db_path} was created successfully but is empty. No rbh-hits were found.")
        return

    m8_outfile_raw = os.path.join(temp_dir, f'{genome1}_vs_{genome2}.m8.raw')
    convert_command = f'mmseqs convertalis {db_1_path} {db_2_path} {result_db_path} {m8_outfile_raw} ' \
                      f'--format-output {consts.MMSEQS_OUTPUT_FORMAT} --search-type 1 --threads 1 -v 1'
    logger.info(f'Calling: {convert_command}')
    subprocess.run(convert_command, shell=True, check=True)

    logger.info(f"{m8_outfile_raw} was created successfully. Adding 'score' column to it...")
    # Add 'score' column to mmseqs output
    rbh_df = pd.read_csv(m8_outfile_raw, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(rbh_df)
    rbh_df = rbh_df[['query', 'target', 'score']]

    if use_parquet:
        rbh_df.to_parquet(output_rbh_path, index=False)
    else:
        rbh_df.to_csv(output_rbh_path, index=False)
    logger.info(f"{output_rbh_path} was created successfully.")

    # Calculate statistics of scores
    scores_statistics = {'mean': statistics.mean(rbh_df['score']), 'sum': sum(rbh_df['score']),
                         'number of records': len(rbh_df['score'])}
    with open(output_statistics_path, 'w') as fp:
        json.dump(scores_statistics, fp)

    logger.info(f"{output_statistics_path} was created successfully.")

    # Calculate max rbh score for each gene
    genome1_max_scores = rbh_df.groupby(['query']).max(numeric_only=True)['score']
    genome2_max_scores = rbh_df.groupby(['target']).max(numeric_only=True)['score']
    genome1_max_scores.to_csv(output_genome1_max_scores, index_label='gene', header=['max_rbh_score'])
    genome2_max_scores.to_csv(output_genome2_max_scores, index_label='gene', header=['max_rbh_score'])

    logger.info(f"{output_genome1_max_scores} and {output_genome2_max_scores} were created successfully.")


def search_rbhs_in_all_pairs(logger, job_input_path, dbs_dir, rbh_hits_dir, scores_statistics_dir,
                             max_rbh_score_per_gene_dir, temp_dir, identity_cutoff, coverage_cutoff, e_value_cutoff,
                             use_parquet):
    with open(job_input_path) as fp:
        for line in fp:
            genome1, genome2 = line.strip().split()
            search_rbh(logger, genome1, genome2, dbs_dir, rbh_hits_dir, scores_statistics_dir,
                       max_rbh_score_per_gene_dir, temp_dir, identity_cutoff, coverage_cutoff, e_value_cutoff,
                       use_parquet)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', help='path to a file that contains the genome pairs to find rbhs')
    parser.add_argument('dbs_dir', help='path to dbs dir')
    parser.add_argument('rbh_hits_dir', help='path to which the results will be written (blast m8 format)')
    parser.add_argument('scores_statistics_dir', help='path to output dir of score statistics')
    parser.add_argument('max_rbh_score_per_gene_dir', help='')
    parser.add_argument('temp_dir', help='path to temp dir')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--use_parquet', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        search_rbhs_in_all_pairs(logger, args.job_input_path, args.dbs_dir, args.rbh_hits_dir,
                                 args.scores_statistics_dir, args.max_rbh_score_per_gene_dir, args.temp_dir,
                                 args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff, args.use_parquet)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
