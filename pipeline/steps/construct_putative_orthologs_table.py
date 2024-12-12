"""
Running example:
script_name.py /Users/Oren/Dropbox/Projects/microbializer/mock_output/all_recip_hits.csv /Users/Oren/Dropbox/Projects/microbializer/mock_output/putative_orthologs_table.tsv
"""

import os
from sys import argv
import argparse
import logging
import sys
from collections import defaultdict
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts


def construct_table(logger, all_reciprocal_hits_path, putative_orthologs_path, single_ogs_dir, delimiter):
    member_gene_to_group_name = {}
    member_gene_to_strain_name_dict = {}
    group_name_to_member_genes = {}
    next_group_id = 0
    with open(all_reciprocal_hits_path) as f:
        for i, line in enumerate(f):
            if i % 1_000_000 == 0:
                logger.info(f'Reciprocal pair number {i}')
            line_tokens = line.rstrip().split(delimiter)
            if 'score' in line:
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

                    group = f'OG_{next_group_id}'
                    next_group_id += 1
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
                    min_group, max_group = sorted([group1, group2])
                    # remove max_group from dictionary
                    max_group_members = group_name_to_member_genes.pop(max_group)
                    # merge groups
                    for gene in max_group_members:
                        member_gene_to_group_name[gene] = min_group
                    group_name_to_member_genes[min_group].extend(max_group_members)

    # strains sorted lexicographically
    sorted_strains = sorted(set(member_gene_to_strain_name_dict.values()))

    header = delimiter.join(['OG_name'] + sorted_strains) + '\n'
    result = ''
    result += header
    # Ignore the group id in group_name_to_member_genes and set a new group id (in the original group id there are
    # many missing ids due to the merging of groups)
    for group_index, (group, group_genes) in enumerate(group_name_to_member_genes.items()):
        # dictionary that holds for each strain the members of the current group
        strain_to_members = defaultdict(list)
        for member in group_genes:
            strain = member_gene_to_strain_name_dict[member]
            strain_to_members[strain].append(member)
        strain_to_members = {strain: ';'.join(members) for strain, members in strain_to_members.items()}
        group_row_str = delimiter.join([f'OG_{group_index}'] + [strain_to_members.get(strain, '') for strain in sorted_strains])
        result += group_row_str + '\n'
        with open(os.path.join(single_ogs_dir, f'OG_{group_index}.txt'), 'w') as f:
            f.write(header)
            f.write(group_row_str)

    with open(putative_orthologs_path, 'w') as f:
        f.write(result)

    output_dir = os.path.dirname(putative_orthologs_path)
    with open(os.path.join(output_dir, 'num_of_putative_sets.txt'), 'w') as f:
        f.write(f'{len(group_name_to_member_genes)}\n')

    if len(group_name_to_member_genes) == 0:
        raise ValueError('No putative orthogroups were found')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('all_reciprocal_hits_path',
                        help='path to a file with all the reciprocal hits files concatenated')
    parser.add_argument('putative_orthologs_path',
                        help='path to an output file in which the putative orthologs table will be written')
    parser.add_argument('single_ogs_dir', help='')
    parser.add_argument('--delimiter', help='delimiter for the input and output files', default=consts.CSV_DELIMITER)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        construct_table(logger, args.all_reciprocal_hits_path, args.putative_orthologs_path, args.single_ogs_dir, args.delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
