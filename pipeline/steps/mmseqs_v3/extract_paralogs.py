import os
import sys
from sys import argv
import argparse
import logging
import pandas as pd
import traceback
import json
import statistics
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
from auxiliaries.logic_auxiliaries import add_score_column_to_mmseqs_output
from auxiliaries import consts


def extract_paralogs_of_genome(logger, m8_path, genome_name, max_scores_parts_dir, paralogs_dir,
                               max_rbh_scores_unified_dir, scores_statistics_dir, temp_dir, use_parquet, verbose):
    genome_max_rbh_scores_path = os.path.join(max_rbh_scores_unified_dir, f'{genome_name}.csv')
    output_paralogs_filtered_path = os.path.join(paralogs_dir, f'{genome_name}_vs_{genome_name}.m8_filtered')
    score_stats_file = os.path.join(scores_statistics_dir, f'{genome_name}_vs_{genome_name}.stats')

    if os.path.exists(genome_max_rbh_scores_path) and os.path.exists(output_paralogs_filtered_path) and os.path.exists(score_stats_file):
        return

    # Unify all max_rbh_scores files of the genome to one file
    max_scores_files = [f for f in os.listdir(max_scores_parts_dir) if f.split('_max_scores_with_')[0] == genome_name]
    max_scores_dfs = []
    for max_scores_file in max_scores_files:
        max_scores_path = os.path.join(max_scores_parts_dir, max_scores_file)
        max_scores_df = pd.read_csv(max_scores_path, index_col='gene')
        max_scores_dfs.append(max_scores_df)

    max_scores_combined_df = pd.concat(max_scores_dfs)
    max_score_per_gene = max_scores_combined_df.groupby('gene')['max_rbh_score'].max()
    max_score_per_gene.to_csv(genome_max_rbh_scores_path)

    # Filter m8_path to include only potential paralogs
    genome_paralogs_path = os.path.join(temp_dir, f'{genome_name}_vs_{genome_name}.m8')
    cmd = f"awk -F'\\t' '$1 ~ /^{genome_name}/ && $2 ~ /^{genome_name}/' {m8_path} > {genome_paralogs_path}"
    subprocess.run(cmd, shell=True)

    m8_df = pd.read_csv(genome_paralogs_path, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    m8_df = m8_df.loc[m8_df['query'] != m8_df['target']]  # Remove self hits
    add_score_column_to_mmseqs_output(m8_df)

    # Keep only hits that have score higher than the max score of both query and target.
    # If only one of the genes was identified as rbh to a gene in another genome (and thus the other one doesn't
    # have max_score value), that's ok.
    # If both genes were not identified as rbh to genes in another genomes, we keep the pair (because we want to
    # keep also paralogs that don't have orthologs in other genomes).
    m8_df['query_max_score'] = m8_df['query'].map(max_score_per_gene)
    m8_df['target_max_score'] = m8_df['target'].map(max_score_per_gene)
    m8_df = m8_df.loc[((m8_df['score'] >= m8_df['query_max_score']) | (m8_df['query_max_score'].isna())) &
                      ((m8_df['score'] >= m8_df['target_max_score']) | (m8_df['target_max_score'].isna()))]

    if m8_df.empty:
        logger.info(f"No paralogs were found for {genome_name} after filtration.")
        return

    # Remove duplicates by sorting query and subject IDs in each pair and taking only unique pairs
    m8_df['pair'] = m8_df.apply(
        lambda x: tuple(sorted((x['query'], x['target']))), axis=1
    )
    m8_df = m8_df.drop_duplicates(subset='pair')

    if use_parquet:
        m8_df[['query', 'target', 'score']].to_parquet(output_paralogs_filtered_path, index=False)
    else:
        m8_df[['query', 'target', 'score']].to_csv(output_paralogs_filtered_path, index=False)
    logger.info(f"{output_paralogs_filtered_path} was created successfully.")

    scores_statistics = {'mean': statistics.mean(m8_df['score']), 'sum': sum(m8_df['score']),
                         'number of records': len(m8_df['score'])}
    with open(score_stats_file, 'w') as fp:
        json.dump(scores_statistics, fp)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('m8_path', help='')
    parser.add_argument('genome_name', help='')
    parser.add_argument('max_scores_parts_dir', help='')
    parser.add_argument('paralogs_dir', help='')
    parser.add_argument('max_rbh_scores_unified_dir', help='')
    parser.add_argument('scores_statistics_dir', help='')
    parser.add_argument('temp_dir', help='')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--use_parquet', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_paralogs_of_genome(logger, args.m8_path, args.genome_name, args.max_scores_parts_dir, args.paralogs_dir,
                                   args.max_rbh_scores_unified_dir, args.scores_statistics_dir, args.temp_dir,
                                   args.use_parquet, args.verbose)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
