import sys
import argparse
import logging

from auxiliaries.pipeline_auxiliaries import load_header2sequences_dict


def filter_out_plasmids(input_genome_path, output_genome_path):
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
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file - which is the input file without plasmids')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    filter_out_plasmids(args.genome_path, args.output_path)
