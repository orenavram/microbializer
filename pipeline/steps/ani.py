from sys import argv
import argparse
import logging
import subprocess
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def reconstruct_msa(logger, query_genome_path, all_genomes_reference_path, output_dir):
    output_path = f'{os.path.join(output_dir, os.path.splitext(os.path.basename(query_genome_path))[0])}_to_all.tsv'
    # No ANI output is reported for a genome pair if ANI value is much below 80% (https://github.com/ParBLiSS/FastANI)
    cmd = f'fastANI -q {query_genome_path} --rl {all_genomes_reference_path} -o {output_path}'
    logger.info(f'Starting fastANI. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('query_genome_path', help='path to a fasta genome file')
    parser.add_argument('all_genomes_reference_path', help='path to a file with a list of all genomes paths')
    parser.add_argument('output_dir', help='path to the output directory')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        reconstruct_msa(logger, args.query_genome_path, args.all_genomes_reference_path, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            f.write(f'Internal Error in {__file__}: {e}\n')
