"""
script_name.py /Users/Oren/Dropbox/Projects/microbializer/mock_output/01_ORFs "GCF_000006945.2_ASM694v2_cds_from_genomic,GCF_000007545.1_ASM754v1_cds_from_genomic,GCF_000008105.1_ASM810v1_cds_from_genomic,GCF_000009505.1_ASM950v1_cds_from_genomic,GCF_000195995.1_ASM19599v1_cds_from_genomic" "lcl_NC_003197.2_cds_NP_459006.1_1_2,lcl_NC_004631.1_cds_WP_001575544.1_1_2,lcl_NC_006905.1_cds_WP_001575544.1_1_2,lcl_NC_011294.1_cds_WP_001575544.1_1_2,lcl_NC_003198.1_cds_NP_454611.1_1_2" /Users/Oren/Dropbox/Projects/microbializer/mock_output/12_orthologs_sets_sequences/og_0.12_orthologs_sets_sequences
"""

import os
from sys import argv
import argparse
import logging
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts


def get_sequence_by_ortholog_name(fasta_path, ortholog_name):
    with open(fasta_path) as f:
        sequence = ''
        flag = False
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if sequence != '':
                    # finished aggregating relevant sequence
                    return sequence
                header = line.lstrip('>')
                if header.startswith(ortholog_name):
                    # gene name was found! start aggregating sequence
                    flag = True
            elif flag:
                # previous header is $ortholog_name! so we are now aggregating its sequence
                sequence += line
        if sequence != '':
            # LAST record was the relevant sequence
            return sequence
    raise ValueError(f'{ortholog_name} does not exist in {fasta_path}!')


def get_orthologs_group_sequences(logger, orfs_dir, strain_name_to_ortholog_name, strains):
    result = ''

    strain_to_strain_orfs_path_dict = {}
    for orfs_file in os.listdir(orfs_dir):
        strain = os.path.splitext(orfs_file)[0]
        strain_to_strain_orfs_path_dict[strain] = os.path.join(orfs_dir, orfs_file)

    for strain in strains:
        ortholog_name = strain_name_to_ortholog_name[strain]
        if ortholog_name != '':
            # current strain has a member in this cluster
            if strain_to_strain_orfs_path_dict.get(strain) is not None:
                orfs_path = strain_to_strain_orfs_path_dict[strain]
                ortholog_sequence = get_sequence_by_ortholog_name(orfs_path, ortholog_name)
                result += f'>{strain}\n{ortholog_sequence}\n'
            else:
                logger.error(f'Could not extract {strain_name_to_ortholog_name[strain]} ortholog of strain {strain} as '
                             f'its ORFs file does not exist at {orfs_dir} (probably ORFs sequence extraction for was '
                             f'failed due to multiple contigs in the corresponding genomic file. Try to grep "failed" '
                             f'on ORFs extraction ER log files)')

    return result


def extract_orfs(logger, sequences_dir, final_table_header_line, cluster_members_line,
                 cluster_name, output_path, delimiter):
    strains = final_table_header_line.rstrip().split(delimiter)
    cluster_members = cluster_members_line.rstrip().split(delimiter)
    strain_name_to_ortholog_name = dict(zip(strains, cluster_members))
    orthologs_sequences = get_orthologs_group_sequences(logger, sequences_dir, strain_name_to_ortholog_name, strains)
    if not orthologs_sequences:
        logger.error(f'Failed to extract any sequence for {cluster_name}.')
        return

    with open(os.path.join(output_path), 'w') as f:
        f.write(orthologs_sequences)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('sequences_dir', help='path to a directory with the bacterial gene sequences (aka ORFs)')
    parser.add_argument('final_table_header', help='string that is the header of the final table')
    parser.add_argument('cluster_members', help='string that is the cluster members that is handled'
                                                '(a row from the final orthologs table)')
    parser.add_argument('cluster_name', help='the name of orthologs group being extracted')
    parser.add_argument('output_path', help='path to an output directory (aka orthologs sets sequences)')
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=consts.CSV_DELIMITER)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_orfs(logger, args.sequences_dir, args.final_table_header, args.cluster_members,
                     args.cluster_name, args.output_path, args.delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
