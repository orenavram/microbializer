import subprocess
import sys
from sys import argv
import argparse
import logging
import os
import traceback
import mmap

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def find_genes(logger, genome, output_dir):
    """
        input:path to fasta file with prokaryotic genome to be analyzed
        output: protein-coding gene prediction for input genome
    """
    fasta_file_prefix = os.path.splitext(os.path.basename(genome))[0]
    orfs_output_file_name = f'{fasta_file_prefix}.fna'
    orfs_output_file_path = os.path.join(output_dir, orfs_output_file_name)
    cmd = f'prodigal -i "{genome}" -d {orfs_output_file_path}'
    logger.info(f'Starting prodigal. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)

    if not os.path.exists(orfs_output_file_path) or os.stat(orfs_output_file_path).st_size == 0:
        raise Exception(f'Could not extract ORFs for {fasta_file_prefix}')
    with open(orfs_output_file_path, 'rb', 0) as orf_f, mmap.mmap(orf_f.fileno(), 0, access=mmap.ACCESS_READ) as s:
        if s.find(b'>') == -1:
            raise Exception(f'{fasta_file_prefix} does not contain any ORFs')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_dir', help='path to output dir')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        find_genes(logger, args.genome_path, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
