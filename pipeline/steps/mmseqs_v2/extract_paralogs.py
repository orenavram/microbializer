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


def extract_paralogs(logger, m8_path, genome_name, max_scores_parts_dir, paralogs_dir,
                     max_rbh_scores_unified_dir, scores_statistics_dir, error_file_path):
    genome_max_rbh_scores_path = os.path.join(max_rbh_scores_unified_dir, f'{genome_name}.csv')
    output_paralogs_raw_path = os.path.join(paralogs_dir, f'{genome_name}_vs_{genome_name}.m8')
    output_paralogs_filtered_path = os.path.join(paralogs_dir, f'{genome_name}_vs_{genome_name}_filtered.m8')
    score_stats_file = os.path.join(scores_statistics_dir, f'{genome_name}_vs_{genome_name}.stats')

    if os.path.exists(genome_max_rbh_scores_path) and os.path.exists(output_paralogs_raw_path) \
            and os.path.exists(output_paralogs_filtered_path) and os.path.exists(score_stats_file):
        return

    m8_df = pd.read_csv(m8_path)

    # Unify all max_rbh_scores files of the genome to one file
    max_score_per_gene = pd.Series(dtype=float)
    max_scores_files = [f for f in os.listdir(max_scores_parts_dir) if f.split('_')[0] == genome_name]
    for max_scores_file in max_scores_files:
        try:
            max_scores_df = pd.read_csv(os.path.join(max_scores_parts_dir, max_scores_file))
            max_score_per_gene = max_score_per_gene.combine(max_scores_df, max_with_nan)
        except Exception as e:
            logger.exception(f'Error while processing {max_scores_file} in method calculate_max_rbh_score_per_gene: {e}')
            fail(logger, f'Error while processing {max_scores_file} in method calculate_max_rbh_score_per_gene: {e}', error_file_path)

    max_score_per_gene.to_csv(genome_max_rbh_scores_path, index_label='gene', header=['max_ortholog_score'])

    # Filter m8_df to include only potential paralogs
    m8_df = m8_df[(m8_df['query_genome'] == genome_name) & (m8_df['target_genome'] == genome_name)]
    m8_df.to_csv(output_paralogs_raw_path, index=False)

    # Keep only hits that have score higher than the max score of both query and target.
    # If only one of the genes was identified as rbh to a gene in another genome (and thus the other one doesn't
    # have max_score value), that's ok.
    # If both genes were not identified as rbh to genes in another genomes, we keep the pair (because we want to
    # keep also paralogs that don't have orthologs in other genomes).
    m8_df['query_max_score'] = m8_df['query'].map(max_score_per_gene)
    m8_df['target_max_score'] = m8_df['target'].map(max_score_per_gene)
    m8_df = m8_df.loc[((m8_df['score'] >= m8_df['query_max_score']) | (m8_df['query_max_score'].isna())) &
                      ((m8_df['score'] >= m8_df['target_max_score']) | (m8_df['target_max_score'].isna()))]

    # Keep only 1 record for each genes pair
    m8_df = m8_df.loc[m8_df['query'] < m8_df['target']]

    if not m8_df.empty:
        m8_df[['query', 'target', 'score']].to_csv(output_paralogs_filtered_path, index=False)
        logger.info(f"{output_paralogs_filtered_path} was created successfully.")

        scores_statistics = {'mean': statistics.mean(m8_df['score']), 'sum': sum(m8_df['score']),
                             'number of records': len(m8_df['score'])}
        with open(score_stats_file, 'w') as fp:
            json.dump(scores_statistics, fp)
    else:
        logger.info(f"No paralogs were found for {genome_name} after filtration.")


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('m8_path', help='')
    parser.add_argument('genome_name', help='path to protein fasta')
    parser.add_argument('max_scores_parts_dir', help='')
    parser.add_argument('paralogs_dir', help='')
    parser.add_argument('max_rbh_scores_unified_dir', help='')
    parser.add_argument('scores_statistics_dir', help='')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_paralogs(logger, args.m8_path, args.genome_name, args.max_scores_parts_dir, args.paralogs_dir,
                         args.max_rbh_scores_unified_dir, args.scores_statistics_dir, args.error_file_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
