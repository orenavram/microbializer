"""
In this module, "core OGs" refer to OGs that contain a gene from each strain, ignoring the user-given argument
core_minimal_percentage. "core genes" are genes that are part of "core OGs".
"""

from pathlib import Path
import sys
import argparse
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step

ORFS_FILE_HEADER_DELIMITER = ' # '
ORFS_FILE_HEADER_GENE_INDEX = 0
ORFS_FILE_HEADER_ORIENTATION_INDEX = 3


def get_genes_info_dicts(fasta_path):
    """
    :param fasta_path: a path to a FASTA file.
           If the file is no in PRODIGAL's output format, default params should be set accordingly.
    :return: two dictionares:
             1. gene_name_to_location: where this gene sits on the genome (first, second, etc..)
             2. gene_name_to_orientation: whether this gene is in sense (1) or anti-sense (-1) orientation
    """
    gene_name_to_location = {}
    gene_name_to_orientation = {}

    location = 0
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            '''
            header example: 
            >Ecoli_APEC_O1_gi_117622295_ref_NC_008563.1__1 # 317 # 2779 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.528
            '''
            line_tokens = line[1:].split(ORFS_FILE_HEADER_DELIMITER)  # returns header without ">" !
            gene_name = line_tokens[ORFS_FILE_HEADER_GENE_INDEX].split(' ')[0]
            orientation = line_tokens[ORFS_FILE_HEADER_ORIENTATION_INDEX]
            gene_name_to_location[gene_name] = location
            gene_name_to_orientation[gene_name] = int(orientation)
            location += 1

    return gene_name_to_location, gene_name_to_orientation


def remove_non_core_genes(genome_name_to_gene_name_to_location, genome_name_to_gene_name_to_orientation,
                          reference_genome_name, genome_names, orthologs_table_path):
    '''
    :param genome_name_to_gene_name_to_location:
    :param genome_name_to_gene_name_to_orientation:
    :param reference_genome_name:
    :param genome_names:
    :param orthologs_table_path:
    :return: a dictionary that maps between a ref gene and its og members and the core genome size
             and changes IN-PLACE genome_name_to_gene_name_to_location & genome_name_to_gene_name_to_orientation
    '''
    ref_gene_to_OG = {}
    genome_name_to_core_genes = {}
    for genome_name in genome_names:
        genome_name_to_core_genes[genome_name] = set()
    with open(orthologs_table_path) as f:
        f.readline()  # skip table's header
        for line in f:
            og = line.rstrip().split(',')[1:]  # skip OG_name
            if not all(og):  # at least one og is missing, i.e., not a core gene
                for genome_name, gene_name in zip(genome_names, og):
                    # discard og members (across participating genomes)
                    genome_name_to_gene_name_to_location[genome_name].pop(gene_name, None)
                    genome_name_to_gene_name_to_orientation[genome_name].pop(gene_name, None)
            else:
                ref_gene = og[0]
                for gene_name, genome_name in zip(og, genome_names):
                    if gene_name in genome_name_to_core_genes[genome_name]:
                        raise ValueError('gene_name in genome_name_to_core_genes[genome_name]!!\n'
                                         'gene_name={gene_name}\n'
                                         'genome_name={genome_name}\n'
                                         'genome_name_to_core_genes[genome_name]={genome_name_to_core_genes[genome_name]}\n')
                    genome_name_to_core_genes[genome_name].add(gene_name)

                if ref_gene in ref_gene_to_OG:
                    raise ValueError('ref_gene in ref_gene_to_OG!!\n'
                                     'ref_gene={ref_gene}\n'
                                     'ref_gene_to_OG={ref_gene_to_OG}\n')
                if ref_gene in genome_name_to_core_genes[reference_genome_name]:
                    ref_gene_to_OG[ref_gene] = og

    # remove "gaps" between locations, i.e., 0-3-4-10 -> 0-1-2-3
    for genome_name, gene_name_to_location in genome_name_to_gene_name_to_location.items():
        location = 0
        for gene_name in sorted(gene_name_to_location, key=gene_name_to_location.get):
            if gene_name not in genome_name_to_core_genes[genome_name]:
                genome_name_to_gene_name_to_location[genome_name].pop(gene_name, None)
                genome_name_to_gene_name_to_orientation[genome_name].pop(gene_name, None)
                continue
            gene_name_to_location[gene_name] = location
            location += 1

    core_genome_size = len(genome_name_to_gene_name_to_location[reference_genome_name])

    return genome_name_to_core_genes, ref_gene_to_OG, core_genome_size


def remove_paralogs(orthologs_table_path, genome_names, no_paralogs_path):
    # Create OG table withhout paralogs

    orthologs_df = pd.read_csv(orthologs_table_path)

    # Remove all rows that have only 1 strain in them (which are orphan genes that have single or multiple copies
    # in the same genome). Those rows will have not-nan values in 2 columns (the OG column and the strain column)
    orthologs_df = orthologs_df[orthologs_df.count(axis=1) > 2]

    # In all other rows, keep max 1 gene per genome
    for index, row in orthologs_df.iterrows():
        for strain in genome_names:
            if pd.notnull(row[strain]) and ';' in row[strain]:
                row[strain] = row[strain].split(';')[0]

    orthologs_df.to_csv(no_paralogs_path, index=False)


def get_genome_numeric_representation(logger, orthologs_table_path, ORFs_dir_path, output_path, tmp_dir):
    with open(orthologs_table_path) as f:
        genome_names = f.readline().rstrip().split(',')[1:]  # skip OG_name

    genome_name_to_gene_name_to_location = {}
    genome_name_to_gene_name_to_orientation = {}
    reference_genome_name = genome_names[0]  # the numbers will be set with respect to this (arbitrary) genome
    for genome_name in genome_names:
        gene_name_to_location, gene_to_orientation = get_genes_info_dicts(ORFs_dir_path / f'{genome_name}.fna')
        genome_name_to_gene_name_to_location[genome_name] = gene_name_to_location
        genome_name_to_gene_name_to_orientation[genome_name] = gene_to_orientation

    no_paralogs_path = tmp_dir / 'no_paralogs.csv'
    remove_paralogs(orthologs_table_path, genome_names, no_paralogs_path)

    genome_name_to_core_genes, ref_gene_to_OG, core_genome_size = \
        remove_non_core_genes(genome_name_to_gene_name_to_location, genome_name_to_gene_name_to_orientation,
                              reference_genome_name, genome_names, no_paralogs_path)

    gene_index = 0
    genome_name_to_numeric_genome = {}
    for genome_name in genome_names:
        genome_name_to_numeric_genome[genome_name] = [None] * core_genome_size

    for ref_gene_name in genome_name_to_gene_name_to_location[reference_genome_name]:
        gene_index += 1

        if ref_gene_name not in genome_name_to_core_genes[reference_genome_name]:
            continue
        og = ref_gene_to_OG[ref_gene_name]
        for gene_name, genome_name in zip(og, genome_names):
            location = genome_name_to_gene_name_to_location[genome_name][gene_name]
            orientation = genome_name_to_gene_name_to_orientation[genome_name][gene_name]
            genome_name_to_numeric_genome[genome_name][location] = str(gene_index * orientation)

    with open(output_path, 'w') as f:
        for genome_name, numeric_genome in genome_name_to_numeric_genome.items():
            numeric_genome = ','.join(numeric_genome)
            f.write(f'>{genome_name}\n{numeric_genome}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', type=Path, help='A path to an ortholog table')
    parser.add_argument('ORFs_dir_path', type=Path, help='A path to a ORF directory')
    parser.add_argument('output_path', type=Path, help='A path to which the numeric core genomes will be written')
    parser.add_argument('tmp_dir', type=Path, help='A path to tmp dir')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, get_genome_numeric_representation, args.orthologs_table_path, args.ORFs_dir_path, args.output_path,
             args.tmp_dir)
