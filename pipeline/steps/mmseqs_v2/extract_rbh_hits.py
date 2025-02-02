import os
import sys
from sys import argv
import argparse
import logging
import pandas as pd
import traceback
import json
import statistics
import dask.dataframe as dd
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, add_default_step_args, str_to_bool


def extract_rbh_hits_of_pair(logger, m8_df, genome1, genome2, rbh_hits_dir, scores_statistics_dir,
                             max_rbh_score_per_gene_dir, temp_dir, use_parquet, verbose):
    output_rbh_path = os.path.join(rbh_hits_dir, f'{genome1}_vs_{genome2}.m8')
    output_statistics_path = os.path.join(scores_statistics_dir, f'{genome1}_vs_{genome2}.stats')
    output_genome1_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome1}_max_scores_with_{genome2}.csv')
    output_genome2_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome2}_max_scores_with_{genome1}.csv')

    if os.path.exists(output_rbh_path) and os.path.exists(output_statistics_path) \
            and os.path.exists(output_genome1_max_scores) and os.path.exists(output_genome2_max_scores):
        return

    genome1_to_2_df = m8_df[(m8_df['query_genome'] == genome1) & (m8_df['target_genome'] == genome2)]
    genome2_to_1_df = m8_df[(m8_df['query_genome'] == genome2) & (m8_df['target_genome'] == genome1)]
    genome1_to_2_df = genome1_to_2_df[['query', 'target', 'score']]
    genome2_to_1_df = genome2_to_1_df[['query', 'target', 'score']]

    if verbose:
        genome1_to_2_df = genome1_to_2_df.sort_values(by=['query', 'target']).reset_index(drop=True)
        genome2_to_1_df = genome2_to_1_df.sort_values(by=['query', 'target']).reset_index(drop=True)
        genome1_to_2_df.to_csv(os.path.join(temp_dir, f'{genome1}_to_{genome2}.m8'), index=False)
        genome2_to_1_df.to_csv(os.path.join(temp_dir, f'{genome2}_to_{genome1}.m8'), index=False)

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

    if reciprocal_best_hits.empty:
        logger.warning(f"No rbh hits were found for {genome1} and {genome2}.")
        return

    # Step 4: Remove duplicates by sorting query and subject IDs in each pair and taking only unique pairs
    reciprocal_best_hits['pair'] = reciprocal_best_hits.apply(
        lambda x: tuple(sorted((x['query_x'], x['target_x']))), axis=1
    )
    unique_rbh = reciprocal_best_hits.drop_duplicates(subset='pair')

    # Step 5: Calculate the average score for each reciprocal best hit
    unique_rbh['average_score'] = int((unique_rbh['score_x'] + unique_rbh['score_y']) / 2)

    # Step 6: Select only relevant columns for final output
    rbh_pairs = unique_rbh[['query_x', 'target_x', 'average_score']]
    rbh_pairs.columns = ['query', 'target', 'score']
    rbh_pairs = rbh_pairs.sort_values(by=['query', 'target']).reset_index(drop=True)

    if use_parquet:
        rbh_pairs.to_parquet(output_rbh_path, index=False)
    else:
        rbh_pairs.to_csv(output_rbh_path, index=False)
    logger.info(f"{output_rbh_path} was created successfully.")

    # Step 7: Calculate statistics of scores
    scores_statistics = {'mean': statistics.mean(rbh_pairs['score']), 'sum': sum(rbh_pairs['score']),
                         'number of records': len(rbh_pairs['score'])}
    with open(output_statistics_path, 'w') as fp:
        json.dump(scores_statistics, fp)
    logger.info(f"{output_statistics_path} was created successfully.")

    # Step 8: Calculate max rbh score for each gene
    genome1_max_scores = rbh_pairs.groupby(['query']).max(numeric_only=True)['score'].reset_index()
    genome1_max_scores.rename(columns={'query': 'gene', 'score': 'max_rbh_score'}, inplace=True)
    genome2_max_scores = rbh_pairs.groupby(['target']).max(numeric_only=True)['score'].reset_index()
    genome2_max_scores.rename(columns={'target': 'gene', 'score': 'max_rbh_score'}, inplace=True)

    if use_parquet:
        genome1_max_scores.to_parquet(output_genome1_max_scores, index=False)
        genome2_max_scores.to_parquet(output_genome2_max_scores, index=False)
    else:
        genome1_max_scores.to_csv(output_genome1_max_scores, index=False)
        genome2_max_scores.to_csv(output_genome2_max_scores, index=False)

    logger.info(f"{output_genome1_max_scores} and {output_genome2_max_scores} were created successfully.")


def extract_rbh_hits(logger, m8_path, rbh_input_path, rbh_hits_dir, scores_statistics_dir, max_rbh_score_per_gene_dir,
                     use_parquet, verbose):
    with open(rbh_input_path, 'r') as f:
        genome_pairs = f.readlines()
        genome_pairs = [pair.strip().split() for pair in genome_pairs]

    temp_dir = os.path.join(rbh_hits_dir, 'tmp')
    os.makedirs(temp_dir, exist_ok=True)

    m8_df = dd.read_parquet(m8_path).compute()
    for genome1, genome2 in genome_pairs:
        extract_rbh_hits_of_pair(logger, m8_df, genome1, genome2, rbh_hits_dir, scores_statistics_dir,
                                 max_rbh_score_per_gene_dir, temp_dir, use_parquet, verbose)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('m8_path', help='a')
    parser.add_argument('rbh_input_path', help='')
    parser.add_argument('rbh_hits_dir', help='')
    parser.add_argument('scores_statistics_dir', help='')
    parser.add_argument('max_rbh_score_per_gene_dir', help='')
    parser.add_argument('--use_parquet', type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        extract_rbh_hits(logger, args.m8_path, args.rbh_input_path, args.rbh_hits_dir,
                         args.scores_statistics_dir, args.max_rbh_score_per_gene_dir, args.use_parquet, args.verbose)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
