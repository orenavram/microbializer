import os
from sys import argv
import argparse
import logging
from Bio import SeqIO
import sys
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def fna_to_faa(logger, nucleotide_path, protein_path, ):
    with open(nucleotide_path, "r") as in_handle, open(protein_path, "w") as out_handle:
        # Iterate through each sequence record in the input file
        for record in SeqIO.parse(in_handle, "fasta"):
            # Translate the DNA sequence into a protein sequence
            translated_record = record.translate(id=True, name=True, description=True)

            # Write the translated record to the output file
            SeqIO.write(translated_record, out_handle, "fasta")

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
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        fna_to_faa(logger, args.nucleotide_path, args.protein_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            f.write(f'Internal Error in {__file__}: {e}\n')
