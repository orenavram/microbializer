import sys
from sys import argv
import argparse
import logging
import os
import traceback
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args


def filter_out_plasmids(logger, input_genome_path, output_genome_path, drop_plasmids, fix_frames):
    """
        input_genome_path: path to an input fasta file with prokaryotic genome
        output_genome_path: path to a filtered genome without plasmids
    """
    if drop_plasmids:
        logger.info(f'Removing plasmids from {input_genome_path}...')
    if fix_frames:
        logger.info(f'Fixing frames of {input_genome_path}...')

    new_records = []
    for record in SeqIO.parse(input_genome_path, 'fasta'):
        if drop_plasmids and ('plasmid' in record.id.lower() or 'plasmid' in record.description.lower()):
            logger.info(f'Dropped plasmid sequence {record.id}')
        else:
            new_records.append(record)

        if fix_frames:
            if 'frame=2' in record.description:
                record.seq = record.seq[1:]
                logger.info(f'Fixed frame of {record.id} (frame=2)')
            elif 'frame=3' in record.description:
                record.seq = record.seq[2:]
                logger.info(f'Fixed frame of {record.id} (frame=3)')

    if not new_records:
        logger.info(f'No records left for {input_genome_path} (probably contained only plasmids)')
        return

    SeqIO.write(new_records, output_genome_path, 'fasta')


def filter_out_plasmids_of_all_files(logger, job_input_path, output_dir, drop_plasmids, fix_frames):
    with open(job_input_path, 'r') as f:
        for line in f:
            genome_path = line.strip()
            genome_file_name = os.path.basename(genome_path)
            output_path = os.path.join(output_dir, genome_file_name)
            filter_out_plasmids(logger, genome_path, output_path, drop_plasmids, fix_frames)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', help='path to a file that contains the genome names to drop plasmids from')
    parser.add_argument('output_dir', help='path to output dir')
    parser.add_argument('--drop_plasmids', action='store_true', help='Drop plasmids from the genome file')
    parser.add_argument('--fix_frames', action='store_true', help='Fix frames of the genome file')
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        filter_out_plasmids_of_all_files(logger, args.job_input_path, args.output_dir, args.drop_plasmids, args.fix_frames)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
