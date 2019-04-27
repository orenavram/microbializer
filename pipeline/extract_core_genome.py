"""
script_name.py and then
/Users/Oren/Downloads/14_aligned_aa_orthologs_groups 21 /Users/Oren/Downloads/11_final_table/strains_names.txt /Users/Oren/Downloads/15_aligned_core_genome/aligned_core_genome.fasta /Users/Oren/Downloads/15_aligned_core_genome/core_ortholog_groups_names.txt --core_minimal_percentage 85
"""

import os

def load_orthologs_group_to_dict(fasta_path):
    header_to_sequence_dict = {}
    with open(fasta_path) as f:
        header = f.readline().lstrip('>').rstrip()
        sequence = ''
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                header_to_sequence_dict[header] = sequence
                header = line.lstrip('>')
                sequence = ''
            else:
                sequence += line

        # don't forget last record!!
        if sequence != '':
            header_to_sequence_dict[header] = sequence

    return header_to_sequence_dict


def is_core_gene(strain_to_gene_dict, num_of_strains, core_minimal_percentage):
    return len(strain_to_gene_dict)/num_of_strains >= core_minimal_percentage/100


def update_core_genome(strain_to_core_genome_dict, strain_to_gene_dict):
    # for strain in strain_to_core_genome_dict:
    #     core_genome_length = len(strain_to_core_genome_dict[strain])
    #     break # all core genomes has the same length
    # else:
    #     # updating core genome for the first time (strain_to_core_genome_dict is empty)
    #     core_genome_length = 0

    for strain in strain_to_gene_dict:
        current_gene_length = len(strain_to_gene_dict[strain])
        break # all alignment genes has the same length

    # strains_to_be_updated = set(strain_to_core_genome_dict).union(set(strain_to_gene_dict))
    for strain in strain_to_core_genome_dict:
        # if strain in strain_to_core_genome_dict:
        #     logger.debug(f'Strain {strain} was already seen.')
        if strain in strain_to_gene_dict:
            logger.debug(f'{strain} has homolog in current og. Updating its core genome')
            strain_to_core_genome_dict[strain] += strain_to_gene_dict[strain]
        else:
            logger.debug(f'{strain} has no homolog in current og. Elongating its core genome by "-"')
            strain_to_core_genome_dict[strain] += '-' * current_gene_length
        # else:
        #     logger.fatal(f'Strain {strain} was encountered for the first time.')
        #     # strain_to_core_genome_dict[strain] = ('-' * core_genome_length) + strain_to_gene_dict[strain]


def extract_core_genome(alignments_path, num_of_strains, core_length_path, strains_names_path, core_genome_path, core_ogs_names_path, core_minimal_percentage):
    with open(strains_names_path) as f:
        strains_names = f.read().rstrip().split('\n')
    strain_to_core_genome_dict = dict.fromkeys(strains_names, '')
    core_ogs = []
    for og_file in os.listdir(alignments_path):  # TODO: consider sorting by og name (currnetly the concatenation is arbitrary)
        strain_to_gene_dict = load_orthologs_group_to_dict(os.path.join(alignments_path, og_file))
        if is_core_gene(strain_to_gene_dict, num_of_strains, core_minimal_percentage):
            logger.info(f'Adding to core genome: {og_file} ({len(strain_to_gene_dict)}/{num_of_strains} >= {core_minimal_percentage}%)')
            update_core_genome(strain_to_core_genome_dict, strain_to_gene_dict)
            core_ogs.append(og_file.split('_')[1])  # e.g., og_2655_aa_mafft.fas
        else:
            logger.info(f'Not a core gene: {og_file} ({len(strain_to_gene_dict)}/{num_of_strains} < {core_minimal_percentage}%)')

    with open(core_genome_path, 'w') as f:
        for strain in sorted(strain_to_core_genome_dict):
            seq = strain_to_core_genome_dict[strain]
            core_genome_length = len(seq)
            f.write(f'>{strain}\n{seq}\n')
            if seq.count('-') == len(seq):
                logger.info("$"*100)
                logger.info(f"seq.count('-')=={seq.count('-')} and len(seq)=={len(seq)}")
                logger.info(f'{strain} has no core genes (core proteome consists only "-" characters)')
                logger.info("$"*100)
                print("$"*100)
                print(f"seq.count('-')=={seq.count('-')} and len(seq)=={len(seq)}")
                print(f'{strain} has no core genes (core proteome consists only "-" characters)')
                print("$"*100)

    with open(core_length_path, 'w') as f:
        f.write(f'{core_genome_length}\n')

    sorted_core_groups_names = sorted(core_ogs, key=int)
    with open(core_ogs_names_path, 'w') as f:
        f.write('\n'.join(f'og_{core_group}' for core_group in sorted_core_groups_names))  # e.g., 2655


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('aa_alignments_path', help='path to a folder where each file is a multiple sequences fasta file')
        parser.add_argument('num_of_strains', help='number of strains in the data', type=int)
        parser.add_argument('strains_names_path', help='path to a file that contains all the strains')
        parser.add_argument('core_genome_path', help='path to an output file in which the core genome will be written')
        parser.add_argument('core_ogs_names_path', help='path to an output file in which the core groups names will be written')
        parser.add_argument('core_length_path', help='path to an output file in which the core genome length')
        parser.add_argument('--core_minimal_percentage', type=float, default=100.0,
                            help='number that represents the required percent that is needed to be considered a core gene. For example: (1) 100 means that for a gene to be considered core, all strains should have a member in the group.\n(2) 50 means that for a gene to be considered core, at least half of the strains should have a member in the group.\n(3) 0 means that every gene should be considered as a core gene.')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        extract_core_genome(args.aa_alignments_path, args.num_of_strains, args.core_length_path, args.strains_names_path,
                            args.core_genome_path, args.core_ogs_names_path, args.core_minimal_percentage)
