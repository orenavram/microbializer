import os
import subprocess
import sys
import time
from sys import argv
import argparse
import logging
import shutil
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, add_default_step_args
from auxiliaries import consts


def too_many_trials(logger, cmd, error_file_path):
    msg = f'Failed to fetch {cmd} command. Could be due to heavy load on our web servers. ' \
          'Please try to re-submit your job in a few minutes or contact us for further information.'
    logger.error(f'Writing to error file in {error_file_path}')
    fail(logger, msg, error_file_path)
    raise OSError(msg)


def run_mmseqs(logger, all_proteins_fasta, output_dir, output_path, identity_cutoff, coverage_cutoff,
               e_value_cutoff, number_of_genomes, cpus, error_file_path, sensitivity):
    tmp_dir = os.path.join(output_dir, f'tmp')

    i = 1
    while not os.path.exists(output_path):
        # when the data set is very big some files are not generated because of the heavy load
        # so we need to make sure they will be generated!
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs easy-search {all_proteins_fasta} {all_proteins_fasta} {output_path} {tmp_dir} ' \
              f'--format-output {consts.MMSEQS_OUTPUT_FORMAT} --min-seq-id {identity_cutoff} -c {coverage_cutoff} ' \
              f'--cov-mode 0 -e {e_value_cutoff} --threads {cpus} -v 1 --comp-bias-corr 0 ' \
              f'--max-seqs {number_of_genomes * 100} -s {sensitivity}'
        logger.info(f'Iteration #{i} - Calling:\n{cmd}')
        subprocess.run(cmd, shell=True, check=True)
        i += 1
        if i == 1000:
            too_many_trials(logger, 'mmseqs easy-rbh', error_file_path)
        time.sleep(1)

        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            if not os.path.exists(output_path):
                tmp_dir = f"{tmp_dir}_try_{i}"

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
    parser.add_argument('--sensitivity', type=float)
    parser.add_argument('--number_of_genomes', type=int, help='Number of genomes in analysis')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use')
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        run_mmseqs(logger, args.all_proteins_fasta, args.output_dir, args.output_path, args.identity_cutoff,
                   args.coverage_cutoff, args.e_value_cutoff, args.number_of_genomes, args.cpus, args.error_file_path,
                   args.sensitivity)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
