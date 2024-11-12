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
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
from auxiliaries.logic_auxiliaries import add_score_column_to_mmseqs_output
from auxiliaries import consts


def extract_rbh_hits(logger, m8_df, genome1, genome2, output_dir, scores_statistics_dir):
    genome1_to_2_df = m8_df[(m8_df['query_genome'] == genome1) & (m8_df['target_genome'] == genome2)]
    genome2_to_1_df = m8_df[(m8_df['query_genome'] == genome2) & (m8_df['target_genome'] == genome1)]

    # Step 1: Identify the best hits from genome1 to genome2
    genome1_to_2_best_hits_df = genome1_to_2_df[genome1_to_2_df['score'] == genome1_to_2_df.groupby('query')['score'].transform('max')]

    # Step 2: Identify the best hits from genome2 to genome1
    genome2_to_1_best_hits_df = genome2_to_1_df[genome2_to_1_df['score'] == genome2_to_1_df.groupby('query')['score'].transform('max')]

    # Step 3: Perform an inner merge on both dataframes to find reciprocal matches
    reciprocal_best_hits = pd.merge(
        genome1_to_2_best_hits_df,
        genome2_to_1_best_hits_df,
        left_on=['query', 'target'],
        right_on=['target', 'query']
    )

    # Step 4: Remove duplicates by sorting query and subject IDs in each pair and taking only unique pairs
    reciprocal_best_hits['pair'] = reciprocal_best_hits.apply(
        lambda x: tuple(sorted((x['query_x'], x['target_x']))), axis=1
    )
    unique_rbh = reciprocal_best_hits.drop_duplicates(subset='pair')

    # Step 5: Calculate the average score for each reciprocal best hit
    unique_rbh['average_score'] = (unique_rbh['score_x'] + unique_rbh['score_y']) / 2

    # Step 6: Select only relevant columns for final output
    rbh_pairs = unique_rbh[['query_x', 'target_x', 'average_score']]
    rbh_pairs.columns = ['query', 'target', 'score']

    rbh_pairs.to_csv(os.path.join(output_dir, f'{genome1}_vs_{genome2}.m8'), sep='\t', index=False)

    score_stats_file = os.path.join(scores_statistics_dir, f'{genome1}_vs_{genome2}.stats')
    scores_statistics = {'mean': statistics.mean(rbh_pairs['score']), 'sum': sum(rbh_pairs['score']),
                         'number of records': len(rbh_pairs['score'])}
    with open(score_stats_file, 'w') as fp:
        json.dump(scores_statistics, fp)


def too_many_trials(logger, cmd, error_file_path):
    msg = f'Failed to fetch {cmd} command. Could be due to heavy load on our web servers. ' \
          'Please try to re-submit your job in a few minutes or contact us for further information.'
    logger.error(f'Writing to error file in {error_file_path}')
    fail(logger, msg, error_file_path)
    raise OSError(msg)


def mmseqs_a_b_c(logger, all_proteins_fasta, strains_names_path, output_dir, scores_statistics_dir, error_file_path,
                      identity_cutoff, coverage_cutoff, e_value_cutoff, cpus):

    # Step a - run mmseqs easy-search on all proteins vs all proteins
    step_a_output_dir = os.path.join(output_dir, '')
    os.makedirs(step_a_output_dir, exist_ok=True)

    tmp_dir = os.path.join(os.path.dirname(step_a_output_dir), f'tmp')
    m8_outfile_raw = os.path.join(step_a_output_dir, 'all_vs_all_raw.m8')
    i = 1
    while not os.path.exists(m8_outfile_raw):
        # when the data set is very big some files are not generated because of the heavy load
        # so we need to make sure they will be generated!
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs easy-rbh {all_proteins_fasta} {all_proteins_fasta} {m8_outfile_raw} {tmp_dir} ' \
              f'--format-output {consts.MMSEQS_OUTPUT_FORMAT} --min-seq-id {identity_cutoff} -c {coverage_cutoff} ' \
              f'--cov-mode 0 -e {e_value_cutoff} --threads {cpus} -v 1'
        logger.info(f'Iteration #{i} - Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)
        i += 1
        if i == 1000:
            too_many_trials(logger, 'mmseqs easy-rbh', error_file_path)
        time.sleep(1)

        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            if not os.path.exists(m8_outfile_raw):
                tmp_dir = f"{tmp_dir}_try_{i}"

    logger.info(f"{m8_outfile_raw} was created successfully. Adding 'score' column to it...")

    # Add 'score' column to mmseqs output
    df = pd.read_csv(m8_outfile_raw, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(df)
    df['query_genome'] = df['query'].str.split(':').str[0]
    df['target_genome'] = df['target'].str.split(':').str[0]

    m8_outfile = os.path.join(step_a_output_dir, 'all_vs_all.m8')
    df.to_csv(m8_outfile, index=False)
    logger.info(f"{m8_outfile} was created successfully.")

    # Step b - extract rbh hits
    df = df[['query', 'query_genome', 'target', 'target_genome', 'score']]




if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('all_proteins_fasta', help='path to a protein fasta')
    parser.add_argument('strains_names_path', help='path to strains names file')
    parser.add_argument('output_dir', help='path to which the results will be written')
    parser.add_argument('scores_statistics_dir', help='path to output dir of score statistics')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('--cpus', type=int)
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        search_all_vs_all(logger, args.all_proteins_fasta, args.strains_names_path, args.output_dir, args.scores_statistics_dir,
                          args.error_file_path, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff, args.cpus)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
