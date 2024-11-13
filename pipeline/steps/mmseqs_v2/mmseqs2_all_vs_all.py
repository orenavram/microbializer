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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
from auxiliaries.logic_auxiliaries import add_score_column_to_mmseqs_output
from auxiliaries import consts


def too_many_trials(logger, cmd, error_file_path):
    msg = f'Failed to fetch {cmd} command. Could be due to heavy load on our web servers. ' \
          'Please try to re-submit your job in a few minutes or contact us for further information.'
    logger.error(f'Writing to error file in {error_file_path}')
    fail(logger, msg, error_file_path)
    raise OSError(msg)


def run_mmseqs(logger, all_proteins_fasta, output_dir, output_path, identity_cutoff, coverage_cutoff,
               e_value_cutoff, cpus, error_file_path):
    tmp_dir = os.path.join(os.path.dirname(output_dir), f'tmp')
    m8_outfile_raw = os.path.join(output_dir, 'all_vs_all_raw.m8')

    i = 1
    while not os.path.exists(m8_outfile_raw):
        # when the data set is very big some files are not generated because of the heavy load
        # so we need to make sure they will be generated!
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs easy-search {all_proteins_fasta} {all_proteins_fasta} {m8_outfile_raw} {tmp_dir} ' \
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

    m8_df = pd.read_csv(m8_outfile_raw, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(m8_df)
    m8_df['query_genome'] = m8_df['query'].str.split(':').str[0]
    m8_df['target_genome'] = m8_df['target'].str.split(':').str[0]

    m8_outfile = os.path.join(output_dir, 'all_vs_all.m8')
    m8_df.to_csv(m8_outfile, index=False)
    logger.info(f"{m8_outfile} was created successfully.")

    m8_df = m8_df[['query', 'query_genome', 'target', 'target_genome', 'score']]
    m8_df.to_csv(output_path, index=False)
    logger.info(f"{output_path} was created successfully.")


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('all_proteins_fasta', help='path to a protein fasta')
    parser.add_argument('output_dir', help='')
    parser.add_argument('output_path', help='')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        run_mmseqs(logger, args.all_proteins_fasta, args.output_dir, args.output_path, args.identity_cutoff,
                   args.coverage_cutoff, args.e_value_cutoff, args.cpus, args.error_file_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
