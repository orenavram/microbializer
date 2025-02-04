import sys
from sys import argv
import argparse
from pathlib import Path
import pandas as pd
import traceback
import json
import statistics
import subprocess
import shutil

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, add_default_step_args, str_to_bool
from auxiliaries.logic_auxiliaries import add_score_column_to_mmseqs_output
from auxiliaries import consts


def search_paralogs(logger, genome_name, dbs_dir, max_scores_parts_dir, paralogs_dir, max_rbh_scores_unified_dir,
                    scores_statistics_dir, temp_dir, identity_cutoff, coverage_cutoff, e_value_cutoff, use_parquet,
                    sensitivity):
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

    genome_db_path = dbs_dir / f'{genome_name}.db'

    tmp_dir = temp_dir / f'tmp_{genome_name}_vs_{genome_name}'
    tmp_dir.mkdir(parents=True, exist_ok=True)

    result_db_path = tmp_dir / f'{genome_name}_vs_{genome_name}.db'
    search_tmp_dir = tmp_dir / f'tmp_search_command'
    search_command = f'mmseqs search {genome_db_path} {genome_db_path} {result_db_path} {search_tmp_dir} ' \
                     f'--min-seq-id {identity_cutoff} -c {coverage_cutoff} --cov-mode 0 -e {e_value_cutoff} --threads 1 ' \
                     f'--search-type 1 --comp-bias-corr 0 -v 1 --alignment-mode 3 -s {sensitivity}'
    logger.info(f'Calling: {search_command}')
    subprocess.run(search_command, shell=True, check=True)

    if result_db_path.stat().st_size == 0:
        logger.info(f"{result_db_path} was created successfully but is empty. No paralogs were found.")
        return

    m8_outfile_raw = tmp_dir / f'{genome_name}_vs_{genome_name}.m8.raw'
    convert_command = f'mmseqs convertalis {genome_db_path} {genome_db_path} {result_db_path} {m8_outfile_raw} ' \
                      f'--format-output {consts.MMSEQS_OUTPUT_FORMAT} --search-type 1 --threads 1 -v 1'
    logger.info(f'Calling: {convert_command}')
    subprocess.run(convert_command, shell=True, check=True)

    logger.info(f"{m8_outfile_raw} was created successfully. Adding 'score' column to it...")
    # Add 'score' column to mmseqs output
    m8_df = pd.read_csv(m8_outfile_raw, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    m8_df = m8_df[m8_df['query'] != m8_df['target']]
    add_score_column_to_mmseqs_output(m8_df)
    m8_df = m8_df[['query', 'target', 'score']]
    m8_df = m8_df.sort_values(by=['query', 'target']).reset_index(drop=True)

    outputs_paralogs_processed = temp_dir / f'{genome_name}_vs_{genome_name}.m8'
    if use_parquet:
        m8_df.to_parquet(outputs_paralogs_processed, index=False)
    else:
        m8_df.to_csv(outputs_paralogs_processed, index=False)
    logger.info(f"{outputs_paralogs_processed} was created successfully.")

    # Delete intermediate files
    shutil.rmtree(tmp_dir, ignore_errors=True)

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


def search_paralogs_in_all_pairs(logger, job_input_path, dbs_dir, max_scores_parts_dir, paralogs_dir,
                                 max_rbh_scores_unified_dir, scores_statistics_dir, temp_dir, identity_cutoff,
                                 coverage_cutoff, e_value_cutoff, use_parquet, sensitivity):
    with open(job_input_path) as fp:
        for line in fp:
            genome_name = line.strip()
            search_paralogs(logger, genome_name, dbs_dir, max_scores_parts_dir, paralogs_dir,
                            max_rbh_scores_unified_dir, scores_statistics_dir, temp_dir, identity_cutoff,
                            coverage_cutoff, e_value_cutoff, use_parquet, sensitivity)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path, help='path to a file that contains the genome names to find paralogs')
    parser.add_argument('dbs_dir', type=Path, help='')
    parser.add_argument('max_scores_parts_dir', type=Path, help='')
    parser.add_argument('paralogs_dir', type=Path, help='')
    parser.add_argument('max_rbh_scores_unified_dir', type=Path, help='')
    parser.add_argument('scores_statistics_dir', type=Path, help='')
    parser.add_argument('temp_dir', type=Path, help='')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('--sensitivity', type=float)
    parser.add_argument('--use_parquet', type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        search_paralogs_in_all_pairs(logger, args.job_input_path, args.dbs_dir, args.max_scores_parts_dir,
                                     args.paralogs_dir, args.max_rbh_scores_unified_dir, args.scores_statistics_dir,
                                     args.temp_dir, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff,
                                     args.use_parquet, args.sensitivity)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
