import sys
import argparse
import logging
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import load_header2sequences_dict, get_job_logger


def filter_out_plasmids(logger, input_genome_path, output_genome_path):
    """
        input_genome_path: path to an input fasta file with prokaryotic genome
        output_genome_path: path to a filtered genome without plasmids
    """
    logger.info(f'Removing plasmids from {input_genome_path}...')
    header2sequences_dict = load_header2sequences_dict(input_genome_path)
    for fasta_header in header2sequences_dict:
        if 'plasmid' in fasta_header:
            logger.info(f'Dropping plasmid sequence {fasta_header}')
            header2sequences_dict.pop(fasta_header)

    if not header2sequences_dict:
        logger.info(f'No records left for {input_genome_path} (probably contained only plasmids)')
        return

    with open(output_genome_path, 'w') as f:
        for fasta_header, sequence in header2sequences_dict.items():
            f.write(f'>{fasta_header}\n{sequence}\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file - which is the input file without plasmids')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    filter_out_plasmids(logger, args.genome_path, args.output_path)
