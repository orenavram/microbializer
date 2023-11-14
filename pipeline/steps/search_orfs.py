import subprocess
import sys
from sys import argv
import argparse
import logging
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def find_genes(logger, genome, output_dir, step_name):
    """
        input:path to fasta file with prokaryotic genome to be analyzed
        output: protein-coding gene prediction for input genome
    """
    fasta_file_prefix = os.path.splitext(genome)[0]
    orfs_output_file_name = f'{fasta_file_prefix}.{step_name}'
    orfs_output_file_path = os.path.join(output_dir, orfs_output_file_name)
    cmd = f'prodigal -i "{genome}" -d {orfs_output_file_path}'
    logger.info(f'Starting prodigal. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_dir', help='path to output dir')
    parser.add_argument('step_name', help='step_name')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        find_genes(logger, args.genome_path, args.output_dir, args.step_name)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
