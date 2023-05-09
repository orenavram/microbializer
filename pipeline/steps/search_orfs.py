import subprocess
import sys
from sys import argv
import argparse
import logging
import os


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def find_genes(logger, genome, output_path):
    """
        input:path to fasta file with prokaryotic genome to be analyzed
        output: protein-coding gene prediction for input genome
    """
    # segments = list(SeqIO.parse(genome, 'fasta'))
    # length = sum(len(segment) for segment in segments)
    # module load prodigal/prodigal-2.6.3
    cmd = f'prodigal -i "{genome}"  -d {output_path}'
    logger.info(f'Starting prodigal. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        find_genes(logger, args.genome_path, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
