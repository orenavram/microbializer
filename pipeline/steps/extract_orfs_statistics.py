from sys import argv
import argparse
import logging
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def extract_orfs_statistics(logger, orf_path, orfs_count_output_path, orfs_gc_output_path):
    num_of_nucleotides = 0
    num_of_GC = 0
    orfs_count = -1
    sequence = ''
    with open(orf_path) as f:
        for line in f:
            if line.startswith('>'):
                orfs_count += 1
                num_of_nucleotides += len(sequence)
                num_of_GC += sequence.count('G') + sequence.count('C')
                sequence = ''
                if orfs_count % 1000 == 0:
                    logger.debug(f'ORFs count is: {orfs_count}')
            else:
                sequence += line.rstrip().upper()

        # don't forget last record!!
        orfs_count += 1
        num_of_nucleotides += len(sequence)
        num_of_GC += sequence.count('G') + sequence.count('C')

        with open(orfs_count_output_path, 'w') as f:
            f.write(f'{orfs_count}\n')  # [f'{count},{orfs_counts[count]}' for count in orfs_counts]))

        with open(orfs_gc_output_path, 'w') as f:
            f.write(f'{num_of_GC / num_of_nucleotides}\n')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('orf_path', help='path to fasta file with orfs')
    parser.add_argument('orfs_count_output_path', help='where to write the number of orfs')
    parser.add_argument('orfs_gc_output_path', help='where to write the gc content')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_orfs_statistics(logger, args.orf_path, args.orfs_count_output_path, args.orfs_gc_output_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
