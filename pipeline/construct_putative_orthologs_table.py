"""
Running example:
script_name.py /Users/Oren/Dropbox/Projects/microbializer/mock_output/all_recip_hits.csv /Users/Oren/Dropbox/Projects/microbializer/mock_output/putative_orthologs_table.tsv
"""

import os


def construct_table(all_reciprocal_hits_path, putative_orthologs_path, delimiter):
    member_gene_to_group_name = {}
    member_gene_to_strain_name_dict = {}
    group_name_to_member_genes = {}
    with open(all_reciprocal_hits_path) as f:
        for line in f:
            line_tokens = line.rstrip().split(delimiter)
            if 'bitscore' in line:
                # new reciprocal hits file starts
                strain1, strain2 = line_tokens[:2]
            else:
                gene1, gene2 = line_tokens[:2]

                gene1_was_seen = gene1 in member_gene_to_strain_name_dict
                gene2_was_seen = gene2 in member_gene_to_strain_name_dict


                if not gene1_was_seen and not gene2_was_seen:
                    # non of the two genes was seen before
                    # update strain name
                    member_gene_to_strain_name_dict[gene1] = strain1
                    member_gene_to_strain_name_dict[gene2] = strain2
                    # take the gene that is associated with the lexicographically lower strain.
                    # Normally, the strains should be different (no paralogs), but if it's the
                    # same strain, take the lexicographically lower gene
                    group = min([gene1, gene2], key=lambda gene: (member_gene_to_strain_name_dict[gene], gene))
                    # add gene1 to group
                    member_gene_to_group_name[gene1] = group
                    # add gene2 to group
                    member_gene_to_group_name[gene2] = group
                    # create group (for the first time)
                    group_name_to_member_genes[group] = [gene1, gene2]

                elif gene1_was_seen and not gene2_was_seen:
                    # update strain name
                    member_gene_to_strain_name_dict[gene2] = strain2
                    # gene1 is already in the table
                    group = member_gene_to_group_name[gene1]
                    # add gene2 to gene1's group
                    member_gene_to_group_name[gene2] = group
                    group_name_to_member_genes[group].append(gene2)

                elif not gene1_was_seen and gene2_was_seen:
                    # update strain name
                    member_gene_to_strain_name_dict[gene1] = strain1

                    # gene2 is already in the table
                    group = member_gene_to_group_name[gene2]
                    # add gene1 to gene2's group
                    member_gene_to_group_name[gene1] = group
                    group_name_to_member_genes[group].append(gene1)

                else:
                    # both genes were already seen; merge groups...
                    group1 = member_gene_to_group_name[gene1]
                    group2 = member_gene_to_group_name[gene2]
                    if group1 == group2:
                        # groups are already merged
                        continue
                    min_group, max_group = sorted([group1, group2], key=lambda gene: (member_gene_to_strain_name_dict[gene], member_gene_to_group_name[gene]))
                    # remove max_group from dictionary
                    max_group_members = group_name_to_member_genes.pop(max_group)
                    # merge groups
                    for gene in max_group_members:
                        member_gene_to_group_name[gene] = min_group
                    group_name_to_member_genes[min_group].extend(max_group_members)

    # strains sorted lexicographically
    sorted_strains = sorted(set(member_gene_to_strain_name_dict.values()))
    putative_orthologs_prefix = os.path.split(putative_orthologs_path)[0]
    with open(os.path.join(putative_orthologs_prefix, 'number_of_strains.txt'), 'w') as f:
        f.write(str(len(sorted_strains)))

    header = delimiter.join(['OG_name']+sorted_strains) + '\n'
    result = ''
    result += header
    # groups sorted by the strain they are associated with (group names are a subset of all the genes..)
    sorted_groups = sorted(group_name_to_member_genes, key=member_gene_to_strain_name_dict.get)
    for group in sorted_groups:
        current_group = group_name_to_member_genes[group]
        # dictionary that holds for each strain the member of the current group if any, o.w., empty str.
        strain_to_member = dict.fromkeys(sorted_strains, '')
        for member in current_group:
            strain = member_gene_to_strain_name_dict[member]
            strain_to_member[strain] = member
        result += delimiter.join([group]+[strain_to_member[strain] for strain in sorted_strains]) + '\n'

    with open(putative_orthologs_path, 'w') as f:
        f.write(result)


if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    logger.info(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('all_reciprocal_hits_path', help='path to a file with all the reciprocal hits files concatenated')
    parser.add_argument('putative_orthologs_path', help='path to an output file in which the putative othologs table will be written')
    parser.add_argument('--delimiter', help='delimiter for the input and output files', default=',')
    args = parser.parse_args()

    construct_table(args.all_reciprocal_hits_path, args.putative_orthologs_path, args.delimiter)