import os
import sys
from sys import argv
import argparse
import logging
import pandas as pd
import traceback
import json
import statistics

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger


def extract_rbh_hits(logger, m8_path, genome1, genome2, rbh_hits_dir, scores_statistics_dir, max_rbh_score_per_gene_dir):
    output_rbh_path = os.path.join(rbh_hits_dir, f'{genome1}_vs_{genome2}.m8')
    output_statistics_path = os.path.join(scores_statistics_dir, f'{genome1}_vs_{genome2}.stats')
    output_genome1_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome1}_max_scores_with_{genome2}.csv')
    output_genome2_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome2}_max_scores_with_{genome1}.csv')

    if os.path.exists(output_rbh_path) and os.path.exists(output_statistics_path) \
            and os.path.exists(output_genome1_max_scores) and os.path.exists(output_genome2_max_scores):
        return

    m8_df = pd.read_csv(m8_path)

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
    rbh_pairs.to_csv(output_rbh_path, sep='\t', index=False)

    # Step 7: Calculate statistics of scores
    scores_statistics = {'mean': statistics.mean(rbh_pairs['score']), 'sum': sum(rbh_pairs['score']),
                         'number of records': len(rbh_pairs['score'])}
    with open(output_statistics_path, 'w') as fp:
        json.dump(scores_statistics, fp)

    # Step 8: Calculate max rbh score for each gene
    genome1_max_scores = rbh_pairs.groupby(['query']).max(numeric_only=True)['score']
    genome1_max_scores.to_csv(output_genome1_max_scores, index_label='gene', header=['max_ortholog_score'])
    genome2_max_scores = rbh_pairs.groupby(['target']).max(numeric_only=True)['score']
    genome2_max_scores.to_csv(output_genome2_max_scores, index_label='gene', header=['max_ortholog_score'])


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('m8_path', help='a')
    parser.add_argument('genome1', help='path to protein fasta')
    parser.add_argument('genome2', help='path to protein fasta')
    parser.add_argument('rbh_hits_dir', help='')
    parser.add_argument('scores_statistics_dir', help='')
    parser.add_argument('max_rbh_score_per_gene_dir', help='')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_rbh_hits(logger, args.m8_path, args.genome1, args.genome2, args.rbh_hits_dir,
                         args.scores_statistics_dir, args.max_rbh_score_per_gene_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)