import sys
from sys import argv
import argparse
import traceback
from pathlib import Path
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args
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
            logger.debug(f'gene {gene_name} was chosen from OG {og_file} to represent strain {strain}')
            strain_to_chosen_gene_sequence_dict[strain] = gene_name_to_sequence_dict[gene_name]

    for strain in strain_to_core_genome_dict:
        if strain in strain_to_chosen_gene_sequence_dict:
            logger.debug(f'{strain} has homolog in og {og_file}. Updating its core genome')
            strain_to_core_genome_dict[strain] += strain_to_chosen_gene_sequence_dict[strain]
        else:
            logger.debug(f'{strain} has no homolog in og {og_file}. Elongating its core genome by "-"')
            strain_to_core_genome_dict[strain] += '-' * og_alignment_length


def extract_core_genome(logger, alignments_path, strains_names_path, core_genome_path, core_length_path,
                        core_minimal_percentage, max_number_of_ogs):
    with open(strains_names_path) as f:
        strains_names = f.read().rstrip().split('\n')
    num_of_strains = len(strains_names)

    output_dir = core_genome_path.parent

    strain_to_core_genome_dict = dict.fromkeys(strains_names, '')
    core_ogs = []
    for og_file in alignments_path.iterdir():  # TODO: consider sorting by og name (currently the concatenation is arbitrary)
        gene_name_to_sequence_dict = {record.id: record.seq for record in SeqIO.parse(og_file, 'fasta')}
        og_alignment_length = len(next(iter(gene_name_to_sequence_dict.values())))
        num_of_strains_in_og = get_num_of_strains_in_og(gene_name_to_sequence_dict)
        if num_of_strains_in_og / num_of_strains >= core_minimal_percentage / 100:  # meaning OG is core
            logger.info(f'Adding to core genome: {og_file.stem} '
                        f'({num_of_strains_in_og}/{num_of_strains} >= {core_minimal_percentage}%)')
            update_core_genome(logger, og_file.stem, gene_name_to_sequence_dict, og_alignment_length, strain_to_core_genome_dict)
            core_ogs.append(og_file.stem.split('_')[1])  # e.g., og_2655 -> 2655
            if max_number_of_ogs and len(core_ogs) == max_number_of_ogs:
                break
        else:
            logger.info(f'Not a core gene: {og_file.stem} '
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
    with open(output_dir / 'core_ortholog_groups_names.txt', 'w') as f:
        f.write('\n'.join(f'OG_{core_group}' for core_group in sorted_core_groups_names))  # e.g., 2655

    with open(output_dir / 'number_of_core_genes.txt', 'w') as f:
        f.write(f'{len(core_ogs)}\n')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('aa_alignments_path', type=Path,
                        help='path to a folder where each file is a multiple sequences fasta file')
    parser.add_argument('strains_names_path', type=Path, help='path to a file that contains all the strains')
    parser.add_argument('core_genome_path', type=Path, help='path to an output file in which the core genome will be written')
    parser.add_argument('core_length_path', type=Path, help='path to an output file in which the core genome length')
    parser.add_argument('--core_minimal_percentage', type=float, default=100.0,
                        help='number that represents the required percent that is needed to be considered a core gene. For example: (1) 100 means that for a gene to be considered core, all strains should have a member in the group.\n(2) 50 means that for a gene to be considered core, at least half of the strains should have a member in the group.\n(3) 0 means that every gene should be considered as a core gene.')
    parser.add_argument('--max_number_of_ogs', type=int, help='maximum number of ogs to add to core genome. None means there is no limit', default=None)
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        extract_core_genome(logger, args.aa_alignments_path, args.strains_names_path, args.core_genome_path,
                            args.core_length_path, args.core_minimal_percentage, args.max_number_of_ogs)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
