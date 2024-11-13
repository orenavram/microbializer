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


def too_many_trials(logger, cmd, error_file_path):
    msg = f'Failed to fetch <i>{cmd}</i> command. Could be due to heavy load on our web servers. ' \
          'Please try to re-submit your job in a few minutes or contact us for further information.'

    # get error_log path
    # e.g., from this aa_db1: /bioseq/data/results/microbializer/159375410340094617808216800611/outputs/02_dbs/SAL_BA5690AA_AS.scaffold_aa
    # into this: /bioseq/data/results/microbializer/159375410340094617808216800611/error.txt
    logger.error(f'Writing to error file in {error_file_path}')
    fail(logger, msg, error_file_path)
    raise OSError(msg)


def search_paralogs(logger, protein_fasta, m8_outfile, genome_max_scores_path, scores_statistics_dir, error_file_path,
                    identity_cutoff, coverage_cutoff, e_value_cutoff):
    """
    input:  protein fasta and a file that contains the max scores of all its genes
    output: mmseqs2 paralogs results file
    """
    strain_name = os.path.splitext(os.path.basename(protein_fasta))[0]
    tmp_dir = os.path.join(os.path.dirname(m8_outfile), f'tmp_{strain_name}')

    i = 1
    while not os.path.exists(m8_outfile):
        # when the data set is very big some files are not generated because of the heavy load
        # so we need to make sure they will be generated!
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs easy-search {protein_fasta} {protein_fasta} {m8_outfile} {tmp_dir} --format-output {consts.MMSEQS_OUTPUT_FORMAT} --min-seq-id {identity_cutoff} -c {coverage_cutoff} --cov-mode 0 -e {e_value_cutoff} --threads 1'
        logger.info(f'Iteration #{i} - Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)
        i += 1
        if i == 1000:
            too_many_trials(logger, 'mmseqs easy-search', error_file_path)
        time.sleep(1)

        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            if not os.path.exists(m8_outfile):
                tmp_dir = f"{tmp_dir}_try_{i}"

    logger.info(f"{m8_outfile} was created successfully. Adding score column and filtering to include only recent paralogs...")
    # Add 'score' column to mmseqs output
    df = pd.read_csv(m8_outfile, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(df)

    # Keep only hits that have score higher than the max score of both query and target.
    # If only one of the genes was identified as homolog to a gene in another genome (and thus the other one doesn't
    # have max_score value), that's ok.
    # If both genes were not identified as homologs to genes in another genome, we keep the pair (because we want to
    # keep also paralogs that don't have orthologs in other genomes).
    genome_max_scores = pd.read_csv(genome_max_scores_path).set_index('gene')['max_ortholog_score']

    df['query_max_score'] = df['query'].map(genome_max_scores)
    df['target_max_score'] = df['target'].map(genome_max_scores)
    df = df.loc[((df['score'] >= df['query_max_score']) | (df['query_max_score'].isna())) &
                ((df['score'] >= df['target_max_score']) | (df['target_max_score'].isna()))]

    # Keep only 1 record for each genes pair
    df = df.loc[df['query'] < df['target']]

    if not df.empty:
        df[['query', 'target', 'score']].to_csv(f'{m8_outfile}_filtered', index=False)
        logger.info(f"{m8_outfile}_filtered was created successfully.")

        score_stats_file = os.path.join(scores_statistics_dir, f'{strain_name}_vs_{strain_name}.stats')
        scores_statistics = {'mean': statistics.mean(df['score']), 'sum': sum(df['score']),
                             'number of records': len(df['score'])}
        with open(score_stats_file, 'w') as fp:
            json.dump(scores_statistics, fp)
    else:
        logger.info(f"No paralogs were found for {strain_name} after filtration.")


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('protein_fasta', help='path to a protein fasta')
    parser.add_argument('output_path', help='path to which the results will be written (blast m8 format)')
    parser.add_argument('genome_max_scores_path',
                        help='path to a file that contains the max ortholog scores of the genes in the fasta')
    parser.add_argument('scores_statistics_dir', help='path to output dir of score statistics')
    parser.add_argument('error_file_path', help='path to which errors are written')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        search_paralogs(logger, args.protein_fasta, args.output_path, args.genome_max_scores_path, args.scores_statistics_dir,
                        args.error_file_path, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
