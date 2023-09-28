from sys import argv
import sys
import argparse
import logging
import os
from Bio import SeqIO
import CodonUsageModified as CodonUsage
import numpy as np
import json
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import get_job_logger

OG_FASTA_HEADER_GENOME_NAME_PATTERN = re.compile(r'.+\((\w+)\)')


def get_genome_to_W_vector(W_dir):
    genome_to_W_vector = {}
    for file_name in os.listdir(W_dir):
        with open(os.path.join(W_dir, file_name), "r") as genome_W_file:
            genome_W = json.load(genome_W_file)
            genome_name = os.path.splitext(file_name)[0]
            genome_to_W_vector[genome_name] = genome_W

    return genome_to_W_vector


def W_vector_to_CAIHandler(W_vector):
    CAI_handler = CodonUsage.CodonAdaptationIndex()
    CAI_handler.set_cai_index(W_vector)
    return CAI_handler


def calculate_cai(OG_dir, W_dir, OG_start_index, OG_stop_index, output_dir, logger):
    genome_to_W_vector = get_genome_to_W_vector(W_dir)
    genome_to_CAIHandler = {genome_name: W_vector_to_CAIHandler(W_vector)
                            for genome_name, W_vector in genome_to_W_vector.items()}
    logger.info(f'Read {len(genome_to_CAIHandler)} W vectors from {W_dir}')

    for OG_index in range(OG_start_index, OG_stop_index + 1):
        cai_info = {}
        OG_path = os.path.join(OG_dir, f'OG_{OG_index}_dna.fas')

        with open(OG_path, 'r') as OG_file:
            for record in SeqIO.parse(OG_file, "fasta"):
                genome_name = OG_FASTA_HEADER_GENOME_NAME_PATTERN.match(record.description).group(1)
                cai_info[genome_name] = genome_to_CAIHandler[genome_name].cai_for_gene(record.seq)

        cai_values = list(cai_info.values())
        cai_info['mean'] = np.mean(cai_values)
        cai_info['std'] = np.std(cai_values)

        output_file_path = os.path.join(output_dir, f'OG_{OG_index}_cai.json')
        with open(output_file_path, 'w') as output_file:
            json.dump(cai_info, output_file)

    logger.info(f'Wrote CAI info for OGs {OG_start_index} - {OG_stop_index} into {output_dir}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('OG_dir', help='path to input Orthologous group directory')
    parser.add_argument('W_dir', help='path to input relative adaptiveness (W vectors) directory')
    parser.add_argument('start', help='starting index of file')
    parser.add_argument('stop', help='stopping index of file (inclusive)')
    parser.add_argument('output_dir', help='path to output dir')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        calculate_cai(args.OG_dir, args.W_dir, int(args.start), int(args.stop), args.output_dir, logger)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
