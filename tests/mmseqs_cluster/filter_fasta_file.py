import sys
from sys import argv
import argparse
import os
from Bio import SeqIO
import traceback

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args


def filter_fasta_file(logger, genes_to_keep_path, input_fasta_path, output_fasta_path):
    with open(genes_to_keep_path, 'r') as fp:
        genes_to_keep = {line.strip() for line in fp}

    with open(input_fasta_path, "r") as infile:
        with open(output_fasta_path, "w") as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                if record.id in genes_to_keep:
                    SeqIO.write(record, outfile, "fasta")


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('genes_to_keep_path', help='path to a file that contains the genes names to keep')
    parser.add_argument('input_fasta_path', help='path to input fasta to filter')
    parser.add_argument('output_fasta_path', help='path to output fasta')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir)

    logger.info(script_run_message)
    try:
        filter_fasta_file(logger, args.genes_to_keep_path, args.input_fasta_path, args.output_fasta_path)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
