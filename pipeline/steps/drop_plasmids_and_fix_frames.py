import sys
from sys import argv
import argparse
import logging
import os
import traceback
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def filter_out_plasmids(logger, input_genome_path, output_genome_path, drop_plasmids, fix_frames):
    """
        input_genome_path: path to an input fasta file with prokaryotic genome
        output_genome_path: path to a filtered genome without plasmids
    """
    logger.info(f'Removing plasmids from {input_genome_path}...')

    new_records = []
    for record in SeqIO.parse(input_genome_path, 'fasta'):
        if drop_plasmids and ('plasmid' in record.id.lower() or 'plasmid' in record.description.lower()):
            logger.info(f'Dropping plasmid sequence {record.id}')
        else:
            new_records.append(record)

        if fix_frames:
            if 'frame=2' in record.description:
                record.seq = record.seq[1:]
            elif 'frame=3' in record.description:
                record.seq = record.seq[2:]

    if not new_records:
        logger.info(f'No records left for {input_genome_path} (probably contained only plasmids)')
        return

    SeqIO.write(new_records, output_genome_path, 'fasta')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file - which is the input file without plasmids')
    parser.add_argument('--drop_plasmids', action='store_true', help='Drop plasmids from the genome file')
    parser.add_argument('--fix_frames', action='store_true', help='Fix frames of the genome file')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        filter_out_plasmids(logger, args.genome_path, args.output_path, args.drop_plasmids, args.fix_frames)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)