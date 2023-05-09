import logging
import datetime
import sys
import os
import argparse


def get_genes_info_dicts(fasta_path, delimiter=' # ', gene_index=0, orientation_index=3):
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
            line_tokens = line[1:].rstrip().split(delimiter)  # returns header without ">" !
            gene_name = line_tokens[gene_index]
            orientation = line_tokens[orientation_index]
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
    ref_gene_to_OGs = {}
    genome_name_to_genes_with_orthologs = {}
    for genome_name in genome_names:
        genome_name_to_genes_with_orthologs[genome_name] = set()
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
                    if gene_name in genome_name_to_genes_with_orthologs[genome_name]:
                        raise ValueError('gene_name in genome_name_to_genes_with_orthologs[genome_name]!!\n'
                                         'gene_name={gene_name}\n'
                                         'genome_name={genome_name}\n'
                                         'genome_name_to_genes_with_orthologs[genome_name]={genome_name_to_genes_with_orthologs[genome_name]}\n')
                    genome_name_to_genes_with_orthologs[genome_name].add(gene_name)

                if ref_gene in ref_gene_to_OGs:
                    raise ValueError('ref_gene in ref_gene_to_OGs!!\n'
                                     'ref_gene={ref_gene}\n'
                                     'ref_gene_to_OGs={ref_gene_to_OGs}\n')
                if ref_gene in genome_name_to_genes_with_orthologs[reference_genome_name]:
                    ref_gene_to_OGs[ref_gene] = og

    # remove "gaps" between locations, i.e., 0-3-4-10 -> 0-1-2-3
    for genome_name in genome_name_to_gene_name_to_location:
        gene_name_to_location = genome_name_to_gene_name_to_location[genome_name]
        location = 0
        for gene_name in sorted(gene_name_to_location, key=gene_name_to_location.get):
            if gene_name not in genome_name_to_genes_with_orthologs[genome_name]:
                genome_name_to_gene_name_to_location[genome_name].pop(gene_name, None)
                genome_name_to_gene_name_to_orientation[genome_name].pop(gene_name, None)
                continue
            gene_name_to_location[gene_name] = location
            location += 1

    core_genome_size = len(genome_name_to_gene_name_to_location[reference_genome_name])

    return genome_name_to_genes_with_orthologs, ref_gene_to_OGs, core_genome_size


def get_genome_numeric_representation(logger, orthologs_table_path, ORFs_dir_path, output_path, output_delimiter=','):
    logger.info(f'{datetime.datetime.now()}: starting to run...')

    with open(orthologs_table_path) as f:
        genome_names = f.readline().rstrip().split(',')[1:]  # skip OG_name

    genome_name_to_gene_name_to_location = {}
    genome_name_to_gene_name_to_orientation = {}
    reference_genome_name = genome_names[0]  # the numbers will be set with respect to this (arbitrary) genome
    for genome_name in genome_names:
        gene_name_to_location, gene_to_orientation = get_genes_info_dicts(
            os.path.join(ORFs_dir_path, f'{genome_name}.02_ORFs'))
        genome_name_to_gene_name_to_location[genome_name] = gene_name_to_location
        genome_name_to_gene_name_to_orientation[genome_name] = gene_to_orientation

    genome_name_to_genes_with_orthologs, ref_gene_to_OGs, core_genome_size = \
        remove_non_core_genes(genome_name_to_gene_name_to_location, genome_name_to_gene_name_to_orientation,
                              reference_genome_name, genome_names, orthologs_table_path)

    gene_index = 0
    genome_name_to_numeric_genome = {}
    for genome_name in genome_names:
        genome_name_to_numeric_genome[genome_name] = [None] * core_genome_size

    for gene_name in genome_name_to_gene_name_to_location[reference_genome_name]:
        gene_index += 1

        if gene_name not in genome_name_to_genes_with_orthologs[reference_genome_name]:
            continue
        og = ref_gene_to_OGs[gene_name]
        for gene_name, genome_name in zip(og, genome_names):
            location = genome_name_to_gene_name_to_location[genome_name][gene_name]
            orientation = genome_name_to_gene_name_to_orientation[genome_name][gene_name]
            genome_name_to_numeric_genome[genome_name][location] = str(gene_index * orientation)

    with open(output_path, 'w') as f:
        for genome_name in genome_name_to_numeric_genome:
            numeric_genome = output_delimiter.join(genome_name_to_numeric_genome[genome_name])
            f.write(f'>{genome_name}\n{numeric_genome}\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', help='A path to an ortholog table (step 11 of microbializer)')
    parser.add_argument('ORFs_dir_path', help='A path to a ORF directory (step 01 of microbializer)')
    parser.add_argument('output_path', help='A path to which the numeric core genomes will be written')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    get_genome_numeric_representation(logger, args.orthologs_table_path, args.ORFs_dir_path.rstrip('/'), args.output_path)
