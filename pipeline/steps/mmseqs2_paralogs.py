import os
import subprocess
import sys
import time
from sys import argv
import argparse
import logging
import shutil
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger
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


def search_paralogs(logger, protein_fasta, m8_outfile, error_file_path, verbosity_level):
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
        logger.info(f'Iteration #{i}: rbh. Result should be at {m8_outfile}')
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs easy-search {protein_fasta} {protein_fasta} {m8_outfile} {tmp_dir} --format-output {consts.MMSEQS_OUTPUT_FORMAT} -v {verbosity_level}'
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)
        i += 1
        if i == 1000:
            too_many_trials(logger, 'mmseqs easy-search', error_file_path)
        time.sleep(1)

    shutil.rmtree(tmp_dir)

    # Remove all rows that contain the same gene as query and subject
    df = pd.read_csv(m8_outfile, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    df = df.loc[df['query'] <= df['target']]
    df.to_csv(m8_outfile, index=False, sep='\t', header=False)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('protein_fasta', help='path to a protein fasta')
    # parser.add_argument('genome_max_scores_path', help='path to a file that contains the max scores of the genes in the fasta')
    parser.add_argument('output_path', help='path to which the results will be written (blast m8 format)')
    parser.add_argument('error_file_path', help='path to which errors are written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        search_paralogs(logger, args.protein_fasta, args.output_path, args.error_file_path, 3 if args.verbose else 1)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')