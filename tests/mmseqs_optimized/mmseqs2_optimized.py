import itertools
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
from multiprocessing import Pool

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
from auxiliaries.logic_auxiliaries import add_score_column_to_mmseqs_output, max_with_nan
from auxiliaries import consts


def extract_paralogs(logger, genome_name, m8_df, max_scores_parts_dir, paralogs_dir,
                     max_rbh_scores_unified_dir, scores_statistics_dir, error_file_path):
    genome_max_rbh_scores_path = os.path.join(max_rbh_scores_unified_dir, f'{genome_name}.csv')
    output_paralogs_raw_path = os.path.join(paralogs_dir, f'{genome_name}_vs_{genome_name}.m8')
    output_paralogs_filtered_path = os.path.join(paralogs_dir, f'{genome_name}_vs_{genome_name}_filtered.m8')
    score_stats_file = os.path.join(scores_statistics_dir, f'{genome_name}_vs_{genome_name}.stats')

    if os.path.exists(genome_max_rbh_scores_path) and os.path.exists(output_paralogs_raw_path) \
            and os.path.exists(output_paralogs_filtered_path) and os.path.exists(score_stats_file):
        return

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


def extract_rbh_hits(logger, m8_df, genome1, genome2, rbh_hits_dir, scores_statistics_dir, max_rbh_score_per_gene_dir):
    output_rbh_path = os.path.join(rbh_hits_dir, f'{genome1}_vs_{genome2}.m8')
    output_statistics_path = os.path.join(scores_statistics_dir, f'{genome1}_vs_{genome2}.stats')
    output_genome1_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome1}_max_scores_with_{genome2}.csv')
    output_genome2_max_scores = os.path.join(max_rbh_score_per_gene_dir, f'{genome2}_max_scores_with_{genome1}.csv')

    if os.path.exists(output_rbh_path) and os.path.exists(output_statistics_path) \
            and os.path.exists(output_genome1_max_scores) and os.path.exists(output_genome2_max_scores):
        return

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


def too_many_trials(logger, cmd, error_file_path):
    msg = f'Failed to fetch {cmd} command. Could be due to heavy load on our web servers. ' \
          'Please try to re-submit your job in a few minutes or contact us for further information.'
    logger.error(f'Writing to error file in {error_file_path}')
    fail(logger, msg, error_file_path)
    raise OSError(msg)


def run_mmseqs_command(logger, all_proteins_fasta, m8_outfile_raw, tmp_dir, error_file_path, identity_cutoff, coverage_cutoff, e_value_cutoff, cpus):
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


def find_orthologs_and_paralogs(logger, all_proteins_fasta, strains_names_path, output_dir, scores_statistics_dir, error_file_path,
                      identity_cutoff, coverage_cutoff, e_value_cutoff, cpus, step_number):
    # Step a - run mmseqs easy-search on all proteins vs all proteins
    step_a_mmseqs_dir = os.path.join(output_dir, f'{step_number}a_mmseqs')
    os.makedirs(step_a_mmseqs_dir, exist_ok=True)

    tmp_dir = os.path.join(os.path.dirname(step_a_mmseqs_dir), f'tmp')
    m8_outfile_raw = os.path.join(step_a_mmseqs_dir, 'all_vs_all_raw.m8')
    run_mmseqs_command(logger, all_proteins_fasta, m8_outfile_raw, tmp_dir, error_file_path, identity_cutoff,
                       coverage_cutoff, e_value_cutoff, cpus)

    m8_df = pd.read_csv(m8_outfile_raw, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(m8_df)
    m8_df['query_genome'] = m8_df['query'].str.split(':').str[0]
    m8_df['target_genome'] = m8_df['target'].str.split(':').str[0]

    m8_outfile = os.path.join(step_a_mmseqs_dir, 'all_vs_all.m8')
    m8_df.to_csv(m8_outfile, index=False)
    logger.info(f"{m8_outfile} was created successfully.")
    m8_outfile_reduced_columns = os.path.join(step_a_mmseqs_dir, 'all_vs_all_reduced_columns.m8')
    m8_df = m8_df[['query', 'query_genome', 'target', 'target_genome', 'score']]
    m8_df.to_csv(m8_outfile_reduced_columns, index=False)
    logger.info(f"{m8_outfile_reduced_columns} was created successfully.")

    # Step b - extract rbh hits
    step_b_rbh_hits_dir = os.path.join(output_dir, f'{step_number}b_rbh_hits')
    os.makedirs(step_b_rbh_hits_dir, exist_ok=True)
    step_b_max_scores_parts_dir = os.path.join(output_dir, f'{step_number}b_max_rbh_scores_parts')
    os.makedirs(step_b_max_scores_parts_dir, exist_ok=True)

    with open(strains_names_path) as f:
        strains_names = f.read().rstrip().split('\n')

    all_rbh_params = []
    for genome1, genome2 in itertools.combinations(strains_names, 2):
        all_rbh_params.append((logger, m8_df, genome1, genome2, step_b_rbh_hits_dir, scores_statistics_dir, step_b_max_scores_parts_dir))

    with Pool(processes=cpus) as pool:
        pool.starmap(extract_rbh_hits, all_rbh_params)

    # Step c - extract paralogs
    step_c_max_rbh_scores_unified_dir = os.path.join(output_dir, f'{step_number}c_max_rbh_scores_unified')
    os.makedirs(step_c_max_rbh_scores_unified_dir, exist_ok=True)
    step_d_paralogs_dir = os.path.join(output_dir, f'{step_number}c_paralogs')
    os.makedirs(step_d_paralogs_dir, exist_ok=True)

    all_paralogs_params = []
    for genome_name in strains_names:
        all_paralogs_params.append((logger, genome_name, m8_df, step_b_max_scores_parts_dir, step_d_paralogs_dir,
                                    step_c_max_rbh_scores_unified_dir, scores_statistics_dir, error_file_path))

    with Pool(processes=cpus) as pool:
        pool.starmap(extract_paralogs, all_paralogs_params)


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
        find_orthologs_and_paralogs(logger, args.all_proteins_fasta, args.strains_names_path, args.output_dir, args.scores_statistics_dir,
                          args.error_file_path, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff, args.cpus)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
