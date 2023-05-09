import sys
from sys import argv
import argparse
import logging
import subprocess
import os
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger


def too_many_trials(logger, cmd, error_file_path):
    msg = f'Failed to fetch <i>{cmd}</i> command. Could be due to heavy load on our web servers. ' \
          'Please try to re-submit your job in a few minutes or contact us for further information.'

    # get error_log path
    # e.g., from this aa_db1: /bioseq/data/results/microbializer/159375410340094617808216800611/outputs/02_dbs/SAL_BA5690AA_AS.scaffold_aa
    # into this: /bioseq/data/results/microbializer/159375410340094617808216800611/error.txt
    logger.error(f'Writing to error file in {error_file_path}')
    fail(logger, msg, error_file_path)
    raise OSError(msg)


def create_mmseq2_DB(logger, fasta_path, output_path, tmp_path, translate, convert2fasta, verbosity_level, cpus=1):
    """
    input:  sequnce to base the DB on
            DB type (nucl/prot)
            path to output file
    output: mmseqs2 DB based on the reference sequence
    """
    # for more details see: mmseqs createdb -h
    i = 1
    while not os.path.exists(tmp_path):
        logger.info(f'Iteration #{i}: createdb. Result should be at {tmp_path}')
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs createdb {fasta_path} {tmp_path} -v {verbosity_level}'  # --dont-split-seq-by-len'
        logger.info(f'Starting mmseqs createdb. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)
        i += 1
        if i == 1000:
            error_file_path = f'{os.path.split(output_path)[0]}/../../error.txt'
            too_many_trials(logger, 'mmseqs createdb', error_file_path)
        time.sleep(1)

    if translate or convert2fasta:
        i = 1
        while not os.path.exists(f'{tmp_path}_aa'):
            logger.info(f'Iteration #{i}: translatenucs. Result should be at {tmp_path}_aa')
            # translate dna db to aa db
            cmd = f'mmseqs translatenucs {tmp_path} {tmp_path}_aa -v {verbosity_level} ' \
                  f'--threads {cpus}'  # default number of threads is 4. If the wrapper did not allocated enough threads the process will crash.
            logger.info(f'Translating dna DB to amino acids. Executed command is:\n{cmd}')
            subprocess.run(cmd, shell=True)
            i += 1
            if i == 1000:
                error_file_path = f'{os.path.split(output_path)[0]}/../../error.txt'
                too_many_trials(logger, 'mmseqs translatenucs', error_file_path)
            time.sleep(1)

    if convert2fasta:  # for dna to aa translation
        i = 1
        while not os.path.exists(f'{output_path}_aa'):
            logger.info(f'Iteration #{i}: convert2fasta. Result should be at {output_path}_aa')
            cmd = f'mmseqs convert2fasta {tmp_path}_aa {output_path}_aa.fas -v {verbosity_level}'  # convert the aa db to a regular fasta
            logger.info(f'Converting amino acids DB to a FASTA format. Executed command is:\n{cmd}')
            subprocess.run(cmd, shell=True)
            i += 1
            if i == 1000:
                error_file_path = f'{os.path.split(output_path)[0]}/../../error.txt'
                too_many_trials(logger, 'mmseqs convert2fasta', error_file_path)
            time.sleep(1)

        intermediate_files = [f'{tmp_path}{suffix}' for suffix in
                              ['', '_aa', '_aa.dbtype', '_aa_h', '_aa_h.dbtype', '_aa_h.index', '_aa.index', '.dbtype',
                               '_h', '_h.dbtype', '_h.index', '.index', '.lookup']]
        # each pair generates 13 intermidiate files!!! lot's of junk once finished
        for file in intermediate_files:
            os.remove(file)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('input_fasta', help='path to input fasta file')
    parser.add_argument('output_prefix', help='path prefix for the DB file(s)')
    parser.add_argument('tmp_prefix', help='path prefix for the tmp file(s)')
    # parser.add_argument('--dbtype', help='database type for creating blast DB', default='nucl')
    parser.add_argument('-t', '--translate', help='whether to translate the dna to aa', action='store_true')
    parser.add_argument('-c', '--convert2fasta', help='whether to convert the dbs to fasta', action='store_true')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        create_mmseq2_DB(logger, args.input_fasta, args.output_prefix, args.tmp_prefix,
                         args.translate, args.convert2fasta, 3 if args.verbose else 1)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
