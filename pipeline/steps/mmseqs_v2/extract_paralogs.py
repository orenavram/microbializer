import sys
from sys import argv
import argparse
from pathlib import Path
import pandas as pd
import traceback
import json
import statistics
import dask.dataframe as dd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, add_default_step_args, str_to_bool


def extract_paralogs_of_genome(logger, m8_df, genome_name, max_scores_parts_dir, paralogs_dir,
                               max_rbh_scores_unified_dir, scores_statistics_dir, temp_dir, use_parquet, verbose):
    genome_max_rbh_scores_path = max_rbh_scores_unified_dir / f'{genome_name}.csv'
    output_paralogs_filtered_path = paralogs_dir / f'{genome_name}_vs_{genome_name}.m8_filtered'
    score_stats_file = scores_statistics_dir / f'{genome_name}_vs_{genome_name}.stats'

    if genome_max_rbh_scores_path.exists() and output_paralogs_filtered_path.exists() and score_stats_file.exists():
        return

    # Unify all max_rbh_scores files of the genome to one file
    max_scores_files = [f for f in max_scores_parts_dir.iterdir() if f.stem.split('_max_scores_with_')[0] == genome_name]
    if max_scores_files:
        max_scores_dfs = []
        for max_scores_file in max_scores_files:
            if use_parquet:
                max_scores_df = pd.read_parquet(max_scores_file)
            else:
                max_scores_df = pd.read_csv(max_scores_file)

            max_scores_dfs.append(max_scores_df)

        max_scores_combined_df = pd.concat(max_scores_dfs)
        max_score_per_gene = max_scores_combined_df.groupby('gene')['max_rbh_score'].max().reset_index()

        if use_parquet:
            max_score_per_gene.to_parquet(genome_max_rbh_scores_path, index=False)
        else:
            max_score_per_gene.to_csv(genome_max_rbh_scores_path, index=False)
        logger.info(f"{genome_max_rbh_scores_path} was created successfully.")

        max_score_per_gene.set_index('gene', inplace=True)
        max_score_per_gene = max_score_per_gene['max_rbh_score']
    else:
        logger.info(f"No max_rbh_scores files were found for {genome_name}.")
        max_score_per_gene = {}

    # Filter m8_df to include only potential paralogs
    m8_df = m8_df[(m8_df['query_genome'] == genome_name) & (m8_df['target_genome'] == genome_name)]
    m8_df = m8_df[['query', 'target', 'score']]
    m8_df = m8_df.sort_values(by=['query', 'target']).reset_index(drop=True)
    output_paralogs_raw_path = temp_dir / f'{genome_name}_vs_{genome_name}.m8'
    if use_parquet:
        m8_df.to_parquet(output_paralogs_raw_path, index=False)
    else:
        m8_df.to_csv(output_paralogs_raw_path, index=False)
    logger.info(f"{output_paralogs_raw_path} was created successfully.")

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
    m8_df = m8_df[['query', 'target', 'score']]
    m8_df = m8_df.sort_values(by=['query', 'target']).reset_index(drop=True)

    if use_parquet:
        m8_df.to_parquet(output_paralogs_filtered_path, index=False)
    else:
        m8_df.to_csv(output_paralogs_filtered_path, index=False)
    logger.info(f"{output_paralogs_filtered_path} was created successfully.")

    scores_statistics = {'mean': statistics.mean(m8_df['score']), 'sum': sum(m8_df['score']),
                         'number of records': len(m8_df['score'])}
    with open(score_stats_file, 'w') as fp:
        json.dump(scores_statistics, fp)

    logger.info(f"{score_stats_file} was created successfully.")


def extract_paralogs(logger, m8_path, genomes_input_path, max_scores_parts_dir, paralogs_dir,
                     max_rbh_scores_unified_dir, scores_statistics_dir, use_parquet, verbose):
    with open(genomes_input_path, 'r') as f:
        genomes = [genome.strip() for genome in f]

    temp_dir = paralogs_dir / 'tmp'
    temp_dir.mkdir(parents=True, exist_ok=True)

    m8_df = dd.read_parquet(m8_path).compute()
    for genome in genomes:
        extract_paralogs_of_genome(logger, m8_df, genome, max_scores_parts_dir, paralogs_dir,
                                   max_rbh_scores_unified_dir, scores_statistics_dir, temp_dir, use_parquet, verbose)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('m8_path', type=Path, help='')
    parser.add_argument('genomes_input_path', type=Path, help='')
    parser.add_argument('max_scores_parts_dir', type=Path, help='')
    parser.add_argument('paralogs_dir', type=Path, help='')
    parser.add_argument('max_rbh_scores_unified_dir', type=Path, help='')
    parser.add_argument('scores_statistics_dir', type=Path, help='')
    parser.add_argument('--use_parquet', type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        extract_paralogs(logger, args.m8_path, args.genomes_input_path, args.max_scores_parts_dir, args.paralogs_dir,
                         args.max_rbh_scores_unified_dir, args.scores_statistics_dir, args.use_parquet, args.verbose)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
