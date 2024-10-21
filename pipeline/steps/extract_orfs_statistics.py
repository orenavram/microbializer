from sys import argv
import argparse
import logging
import os
import sys
import json
from Bio import SeqIO
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def extract_orfs_statistics(logger, orf_path, orfs_statistics_dir):
    total_num_of_nucleotides = 0
    total_num_of_GC = 0
    orfs_count = 0
    for seq_record in SeqIO.parse(orf_path, 'fasta'):
        orfs_count += 1
        total_num_of_nucleotides += len(seq_record)
        total_num_of_GC += seq_record.count('G') + seq_record.count('C') + seq_record.count('g') + seq_record.count('c')

    orfs_statistics = {}
    orfs_statistics['orfs_count'] = orfs_count
    orfs_statistics['gc_content'] = total_num_of_GC / total_num_of_nucleotides

    genome_name = os.path.splitext(os.path.basename(orf_path))[0]
    orfs_statistics_step_name = os.path.basename(orfs_statistics_dir)
    orfs_statistics_file_name = f'{genome_name}.{orfs_statistics_step_name}'
    orfs_statistics_file_path = os.path.join(orfs_statistics_dir, orfs_statistics_file_name)
    with open(orfs_statistics_file_path, 'w') as fp:
        json.dump(orfs_statistics, fp)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('orf_path', help='path to fasta file with orfs')
    parser.add_argument('orfs_statistics_dir', help='output dir of orfs statistics')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_orfs_statistics(logger, args.orf_path, args.orfs_statistics_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
