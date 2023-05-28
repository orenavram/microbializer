import os
import subprocess
import sys
import time
from sys import argv
import argparse
import logging

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


def search_all_vs_all(logger, aa_db1, aa_db2, aln_offsetted_db, tmp_dir, m8_outfile, verbosity_level, cpus=1):
    """
    input:  mmseqs2 DBs
    output: query_vs_reference "mmseqs2 search" results file
    """
    i = 1
    while not os.path.exists(aln_offsetted_db):
        # when the data set is very big some files are not generated because of the heavy load
        # so we need to make sure they will be generated!
        logger.info(f'Iteration #{i}: rbh. Result should be at {aln_offsetted_db}')
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs rbh {aa_db1} {aa_db2} {aln_offsetted_db} {tmp_dir} -v {verbosity_level} --remove-tmp-files ' \
              f'--threads {cpus}'  # default number of threads is 4. If the wrapper did not allocated enough threads the process will crash.
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)
        i += 1
        if i == 1000:
            error_file_path = f'{os.path.split(aa_db1)[0]}/../../error.txt'
            too_many_trials(logger, 'mmseqs rbh', error_file_path)
        time.sleep(1)

    i = 1
    while not os.path.exists(m8_outfile):
        # when the data set is very big some files are not generated because of the heavy load
        logger.info(f'Iteration #{i}: convertalis. Result should be at {m8_outfile}')
        cmd = f'mmseqs convertalis {aa_db1} {aa_db2} {aln_offsetted_db} {m8_outfile} -v {verbosity_level} ' \
              f'--threads {cpus}'  # default number of threads is 4. If the wrapper did not allocated enough threads the process will crash.
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)
        i += 1
        if i == 1000:
            error_file_path = f'{os.path.split(aa_db1)[0]}/../../error.txt'
            too_many_trials(logger, 'mmseqs convertalis', error_file_path)
        time.sleep(1)

    intermediate_files_prefix = os.path.splitext(aln_offsetted_db)[0]
    intermediate_files = [f'{intermediate_files_prefix}{suffix}' for suffix in
                          ['.alnOffsettedDB', '.alnOffsettedDB.dbtype', '.alnOffsettedDB.index']]
    # each pair generates 3 intermidiate files! lot's of junk once finished
    for file in intermediate_files:
        os.remove(file)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('aa_db1', help='path to an aa DB')
    parser.add_argument('aa_db2', help='path to another aa DB')
    parser.add_argument('aln_offsetted_db', help='path to mmseqs2 offsetted alignment DB')
    parser.add_argument('tmp_dir', help='a path to write mmseqs internal files')
    parser.add_argument('output_path', help='path to which the results will be written (blast m8 format)')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        search_all_vs_all(logger, args.aa_db1, args.aa_db2, args.aln_offsetted_db,
                          args.tmp_dir, args.output_path, 3 if args.verbose else 1)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
