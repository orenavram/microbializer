import os
import sys
from sys import argv
import argparse
import logging
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import load_header2sequences_dict, get_job_logger
from auxiliaries.logic_auxiliaries import get_strain_name


def get_num_of_strains_in_og(gene_name_to_sequence_dict):
    strains_in_og = set()
    for gene_name in gene_name_to_sequence_dict.keys():
        strain_name = get_strain_name(gene_name)
        strains_in_og.add(strain_name)

    return len(strains_in_og)


def update_core_genome(logger, og_file, gene_name_to_sequence_dict, og_alignment_length, strain_to_core_genome_dict):
    # gene_name_to_sequence_dict may contain multiple genes for each strain, so here we choose one of them
    # to be included in the core genome.
    strain_to_chosen_gene_sequence_dict = {}
    for gene_name, sequence in gene_name_to_sequence_dict.items():
        strain = get_strain_name(gene_name)
        if strain not in strain_to_chosen_gene_sequence_dict:
            logger.info(f'gene {gene_name} was chosen from OG {og_file} to represent strain {strain}')
            strain_to_chosen_gene_sequence_dict[strain] = gene_name_to_sequence_dict[gene_name]

    for strain in strain_to_core_genome_dict:
        if strain in strain_to_chosen_gene_sequence_dict:
            logger.debug(f'{strain} has homolog in og {og_file}. Updating its core genome')
            strain_to_core_genome_dict[strain] += strain_to_chosen_gene_sequence_dict[strain]
        else:
            logger.debug(f'{strain} has no homolog in og {og_file}. Elongating its core genome by "-"')
            strain_to_core_genome_dict[strain] += '-' * og_alignment_length


def extract_core_genome(logger, alignments_path, num_of_strains, core_length_path, strains_names_path, core_genome_path,
                        core_ogs_names_path, number_of_core_members_path, core_minimal_percentage):
    with open(strains_names_path) as f:
        strains_names = f.read().rstrip().split('\n')
    strain_to_core_genome_dict = dict.fromkeys(strains_names, '')
    core_ogs = []
    for og_file in os.listdir(alignments_path):  # TODO: consider sorting by og name (currently the concatenation is arbitrary)
        gene_name_to_sequence_dict, og_alignment_length = load_header2sequences_dict(
            os.path.join(alignments_path, og_file), get_length=True)
        num_of_strains_in_og = get_num_of_strains_in_og(gene_name_to_sequence_dict)
        if num_of_strains_in_og / num_of_strains >= core_minimal_percentage / 100:  # meaning OG is core
            logger.info(f'Adding to core genome: {og_file} '
                        f'({num_of_strains_in_og}/{num_of_strains} >= {core_minimal_percentage}%)')
            update_core_genome(logger, og_file, gene_name_to_sequence_dict, og_alignment_length, strain_to_core_genome_dict)
            core_ogs.append(og_file.split('_')[1])  # e.g., og_2655_aa_mafft.fas
        else:
            logger.info(f'Not a core gene: {og_file} '
                        f'({num_of_strains_in_og}/{num_of_strains} < {core_minimal_percentage}%)')

    core_genome_length = None
    with open(core_genome_path, 'w') as f:
        for strain in sorted(strain_to_core_genome_dict):
            seq = strain_to_core_genome_dict[strain]
            core_genome_length = len(seq)
            f.write(f'>{strain}\n{seq}\n')
            if seq.count('-') == len(seq):
                logger.error("#" * 100)
                logger.error(f"seq.count('-')=={seq.count('-')} and len(seq)=={len(seq)}")
                logger.error(f'{strain} has no core genes (core proteome consists only "-" characters)')
                logger.error("#" * 100)

    with open(core_length_path, 'w') as f:
        f.write(f'{core_genome_length}\n')

    sorted_core_groups_names = sorted(core_ogs, key=int)
    with open(core_ogs_names_path, 'w') as f:
        f.write('\n'.join(f'OG_{core_group}' for core_group in sorted_core_groups_names))  # e.g., 2655

    with open(number_of_core_members_path, 'w') as f:
        f.write(f'{len(core_ogs)}\n')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('aa_alignments_path',
                        help='path to a folder where each file is a multiple sequences fasta file')
    parser.add_argument('num_of_strains', help='number of strains in the data', type=int)
    parser.add_argument('strains_names_path', help='path to a file that contains all the strains')
    parser.add_argument('core_genome_path', help='path to an output file in which the core genome will be written')
    parser.add_argument('core_ogs_names_path',
                        help='path to an output file in which the core groups names will be written')
    parser.add_argument('core_length_path', help='path to an output file in which the core genome length')
    parser.add_argument('number_of_core_members_path',
                        help='path to an output file in which the number of core genes detected will be written to')
    parser.add_argument('--core_minimal_percentage', type=float, default=100.0,
                        help='number that represents the required percent that is needed to be considered a core gene. For example: (1) 100 means that for a gene to be considered core, all strains should have a member in the group.\n(2) 50 means that for a gene to be considered core, at least half of the strains should have a member in the group.\n(3) 0 means that every gene should be considered as a core gene.')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_core_genome(logger, args.aa_alignments_path, args.num_of_strains, args.core_length_path,
                            args.strains_names_path, args.core_genome_path, args.core_ogs_names_path,
                            args.number_of_core_members_path, args.core_minimal_percentage)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
