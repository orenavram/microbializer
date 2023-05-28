import os
from sys import argv
import argparse
import logging
import Bio.Seq
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def fna_to_faa(logger, nucleotide_path, protein_path):
    with open(nucleotide_path) as f:
        result = ''
        previous_protein = ''
        previous_header = f.readline()
        for line in f:
            if line.startswith('>'):
                result += f'{previous_header}{Bio.Seq.translate(previous_protein)}\n'
                previous_protein = ''
                previous_header = line
            else:
                previous_protein += line.rstrip()
        # don't forget to add last record!
        result += f'{previous_header}{Bio.Seq.translate(previous_protein)}\n'

    with open(protein_path, 'w') as f:
        f.write(result)

    logger.info(f'Translated fatsa file was write successfully to: {protein_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('nucleotide_path',
                        help='A path to a nucleotide fasta file',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('protein_path', help='A path to which the translated dna will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        fna_to_faa(logger, args.nucleotide_path, args.protein_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
