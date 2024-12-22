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
               error_file_path, identity_cutoff, coverage_cutoff, e_value_cutoff, use_parquet):
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

    tmp_dir = os.path.join(rbh_hits_dir, f'tmp_{genome1}_vs_{genome2}')
    m8_outfile_raw = os.path.join(rbh_hits_dir, f'{genome1}_vs_{genome2}.m8.raw')
 
    # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
    db_1_path = os.path.join(dbs_dir, f'{genome1}.db')
    db_2_path = os.path.join(dbs_dir, f'{genome2}.db')
    rbh_command = f'mmseqs rbh {db_1_path} {db_2_path} {{result_db}} {{tmp_dir}} ' \
    f'--min-seq-id {IDENTITY_CUTOFF} ' \
    f'-c {COVERAGE_CUTOFF} --cov-mode 0 -e {E_VALUE_CUTOFF} --threads {{cpus}} --search-type 1 --comp-bias-corr 0'

    cmd = f'mmseqs easy-rbh {protein_fasta_1} {protein_fasta_2} {m8_outfile_raw} {tmp_dir} --format-output {consts.MMSEQS_OUTPUT_FORMAT} --min-seq-id {identity_cutoff} -c {coverage_cutoff} --cov-mode 0 -e {e_value_cutoff} --threads 1 -v 1'
    logger.info(f'Iteration #{i} - Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)


    if os.path.getsize(m8_outfile_raw) == 0:
        logger.info(f"{m8_outfile_raw}  was created successfully but is empty. No rbh-hits were found.")
        return

    logger.info(f"{m8_outfile_raw} was created successfully. Adding 'score' column to it...")
    # Add 'score' column to mmseqs output
    df = pd.read_csv(m8_outfile_raw, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(df)

    if use_parquet:
        df[['query', 'target', 'score']].to_parquet(m8_outfile, index=False)
    else:
        df[['query', 'target', 'score']].to_csv(m8_outfile, index=False)
    logger.info(f"{m8_outfile} was created successfully.")

    scores_statistics = {'mean': statistics.mean(df['score']), 'sum': sum(df['score']),
                         'number of records': len(df['score'])}
    with open(score_stats_file, 'w') as fp:
        json.dump(scores_statistics, fp)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('protein_fasta_1', help='path to a protein fasta')
    parser.add_argument('protein_fasta_2', help='path to another protein fasta')
    parser.add_argument('output_path', help='path to which the results will be written (blast m8 format)')
    parser.add_argument('scores_statistics_dir', help='path to output dir of score statistics')
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
        search_all_vs_all(logger, args.protein_fasta_1, args.protein_fasta_2, args.output_path, args.scores_statistics_dir,
                          args.error_file_path, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff, args.use_parquet)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
