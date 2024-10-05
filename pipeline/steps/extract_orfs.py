"""
script_name.py /Users/Oren/Dropbox/Projects/microbializer/mock_output/01_ORFs "GCF_000006945.2_ASM694v2_cds_from_genomic,GCF_000007545.1_ASM754v1_cds_from_genomic,GCF_000008105.1_ASM810v1_cds_from_genomic,GCF_000009505.1_ASM950v1_cds_from_genomic,GCF_000195995.1_ASM19599v1_cds_from_genomic" "lcl_NC_003197.2_cds_NP_459006.1_1_2,lcl_NC_004631.1_cds_WP_001575544.1_1_2,lcl_NC_006905.1_cds_WP_001575544.1_1_2,lcl_NC_011294.1_cds_WP_001575544.1_1_2,lcl_NC_003198.1_cds_NP_454611.1_1_2" /Users/Oren/Dropbox/Projects/microbializer/mock_output/12_orthologs_sets_sequences/og_0.12_orthologs_sets_sequences
"""

import os
from sys import argv
import argparse
import logging
import sys
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts


def extract_orfs_of_og(logger, strain_to_gene_to_sequence_dict, strain_name_to_genes_names, og_name, output_dir):
    og_sequences = ''
    for strain, genes_names in strain_name_to_genes_names.items():
        if not genes_names:
            continue

        # current strain has members in this cluster
        gene_to_sequence = strain_to_gene_to_sequence_dict[strain]
        og_gene_to_sequence = {gene: gene_to_sequence[gene] for gene in genes_names.split(';')}
        for gene_name, sequence in og_gene_to_sequence.items():
            og_sequences += f'>{gene_name}\n{sequence}\n'

    if not og_sequences:
        logger.error(f'Failed to extract any sequence for {og_name}.')
        return

    output_path = os.path.join(output_dir, f'{og_name}_dna.fas')
    with open(os.path.join(output_path), 'w') as f:
        f.write(og_sequences)


def extract_orfs(logger, orfs_dir, orthogroups_file_path, start_line_index, end_line_index_exclusive, output_dir, delimiter):
    strain_to_gene_to_sequence_dict = {}
    for orfs_file in os.listdir(orfs_dir):
        strain = os.path.splitext(orfs_file)[0]
        gene_to_sequence = {}
        for seq_record in SeqIO.parse(os.path.join(orfs_dir, orfs_file), 'fasta'):
            gene_to_sequence[seq_record.id] = seq_record.seq
        strain_to_gene_to_sequence_dict[strain] = gene_to_sequence

    with open(orthogroups_file_path) as f:
        header_line = f.readline()
        first_delimiter_index = header_line.index(delimiter)
        final_table_header = header_line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"
        strains = final_table_header.rstrip().split(delimiter)

        for i, line in enumerate(f):
            if i < start_line_index or i >= end_line_index_exclusive:
                continue
            first_delimiter_index = line.index(consts.CSV_DELIMITER)
            og_name = line[:first_delimiter_index]
            cluster_members = line.rstrip()[first_delimiter_index + 1:].split(delimiter)  # remove "OG_name" and split
            strain_name_to_genes_names = dict(zip(strains, cluster_members))
            extract_orfs_of_og(logger, strain_to_gene_to_sequence_dict, strain_name_to_genes_names, og_name, output_dir)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('orfs_dir', help='path to a directory with the bacterial gene sequences (aka ORFs)')
    parser.add_argument('orthogroups_file_path', help='path of the orthogroups file')
    parser.add_argument('start_line_index', help='', type=int)
    parser.add_argument('end_line_index_exclusive', help='', type=int)
    parser.add_argument('output_dir', help='path to an output directory (aka orthologs sets sequences)')
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=consts.CSV_DELIMITER)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_orfs(logger, args.orfs_dir, args.orthogroups_file_path, args.start_line_index,
                     args.end_line_index_exclusive, args.output_dir, args.delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
