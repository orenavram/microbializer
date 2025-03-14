import subprocess
import sys
import argparse
from pathlib import Path
import shutil
import pandas as pd
import statistics
import json

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.logic_utils import add_score_column_to_mmseqs_output
from pipeline.auxiliaries.general_utils import str_to_bool
from pipeline.auxiliaries import consts


def search_rbh(logger, genome1, genome2, dbs_dir, rbh_hits_dir, scores_statistics_dir, max_rbh_score_per_gene_dir,
               temp_dir, identity_cutoff, coverage_cutoff, e_value_cutoff, use_parquet, sensitivity):
    """
    input:  2 protein dbs
    output: query_vs_reference "mmseqs2 easy-rbh" results file
    """
    output_rbh_path = rbh_hits_dir / f'{genome1}_vs_{genome2}.m8'
    output_statistics_path = scores_statistics_dir / f'{genome1}_vs_{genome2}.stats'
    output_genome1_max_scores = max_rbh_score_per_gene_dir / f'{genome1}_max_scores_with_{genome2}.csv'
    output_genome2_max_scores = max_rbh_score_per_gene_dir / f'{genome2}_max_scores_with_{genome1}.csv'

    if output_rbh_path.exists() and output_statistics_path.exists() \
            and output_genome1_max_scores.exists() and output_genome2_max_scores.exists():
        return

    tmp_dir = temp_dir / f'tmp_{genome1}_vs_{genome2}'
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
    db_1_path = dbs_dir / f'{genome1}.db'
    db_2_path = dbs_dir / f'{genome2}.db'
    result_db_path = tmp_dir / f'{genome1}_vs_{genome2}.db'
    rbh_command_tmp_dir = tmp_dir / 'tmp_rbh_command'
    rbh_command = f'mmseqs rbh {db_1_path} {db_2_path} {result_db_path} {rbh_command_tmp_dir} --min-seq-id {identity_cutoff} ' \
                  f'-c {coverage_cutoff} --cov-mode 0 -e {e_value_cutoff} --threads 1 --search-type 1 ' \
                  f'--comp-bias-corr 0 -v 1 -s {sensitivity}'
    logger.info(f'Calling: {rbh_command}')
    subprocess.run(rbh_command, shell=True, check=True, capture_output=True, text=True)

    if result_db_path.stat().st_size == 0:
        logger.info(f"{result_db_path} was created successfully but is empty. No rbh-hits were found.")
        return

    m8_outfile_raw = tmp_dir / f'{genome1}_vs_{genome2}.m8.raw'
    convert_command = f'mmseqs convertalis {db_1_path} {db_2_path} {result_db_path} {m8_outfile_raw} ' \
                      f'--format-output {consts.MMSEQS_OUTPUT_FORMAT} --search-type 1 --threads 1 -v 1'
    logger.info(f'Calling: {convert_command}')
    subprocess.run(convert_command, shell=True, check=True, capture_output=True, text=True)

    logger.info(f"{m8_outfile_raw} was created successfully. Adding 'score' column to it...")
    # Add 'score' column to mmseqs output
    rbh_df = pd.read_csv(m8_outfile_raw, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(rbh_df)
    rbh_df = rbh_df[['query', 'target', 'score']]
    rbh_df = rbh_df.sort_values(by=['query', 'target']).reset_index(drop=True)

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
    genome1_max_scores = rbh_df.groupby(['query']).max(numeric_only=True)['score'].reset_index()
    genome1_max_scores.rename(columns={'query': 'gene', 'score': 'max_rbh_score'}, inplace=True)
    genome2_max_scores = rbh_df.groupby(['target']).max(numeric_only=True)['score'].reset_index()
    genome2_max_scores.rename(columns={'target': 'gene', 'score': 'max_rbh_score'}, inplace=True)

    if use_parquet:
        genome1_max_scores.to_parquet(output_genome1_max_scores, index=False)
        genome2_max_scores.to_parquet(output_genome2_max_scores, index=False)
    else:
        genome1_max_scores.to_csv(output_genome1_max_scores, index=False)
        genome2_max_scores.to_csv(output_genome2_max_scores, index=False)

    logger.info(f"{output_genome1_max_scores} and {output_genome2_max_scores} were created successfully.")

    # Delete intermediate files
    shutil.rmtree(tmp_dir, ignore_errors=True)


def search_rbhs_in_all_pairs(logger, job_input_path, dbs_dir, rbh_hits_dir, scores_statistics_dir,
                             max_rbh_score_per_gene_dir, temp_dir, identity_cutoff, coverage_cutoff, e_value_cutoff,
                             use_parquet, sensitivity):
    with open(job_input_path) as fp:
        for line in fp:
            genome1, genome2 = line.strip().split()
            search_rbh(logger, genome1, genome2, dbs_dir, rbh_hits_dir, scores_statistics_dir,
                       max_rbh_score_per_gene_dir, temp_dir, identity_cutoff, coverage_cutoff, e_value_cutoff,
                       use_parquet, sensitivity)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path, help='path to a file that contains the genome pairs to find rbhs')
    parser.add_argument('dbs_dir', type=Path, help='path to dbs dir')
    parser.add_argument('rbh_hits_dir', type=Path, help='path to which the results will be written (blast m8 format)')
    parser.add_argument('scores_statistics_dir', type=Path, help='path to output dir of score statistics')
    parser.add_argument('max_rbh_score_per_gene_dir', type=Path, help='')
    parser.add_argument('temp_dir', type=Path, help='path to temp dir')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('--sensitivity', type=float)
    parser.add_argument('--use_parquet', type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, search_rbhs_in_all_pairs, args.job_input_path, args.dbs_dir, args.rbh_hits_dir,
             args.scores_statistics_dir, args.max_rbh_score_per_gene_dir, args.temp_dir, args.identity_cutoff,
             args.coverage_cutoff, args.e_value_cutoff, args.use_parquet, args.sensitivity)
