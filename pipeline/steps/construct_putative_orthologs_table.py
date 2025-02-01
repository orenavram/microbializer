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
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args
from auxiliaries import consts


def cluster_genes_to_connected_components(logger, normalized_hits_dir, putative_orthologs_path, delimiter):
    member_gene_to_group_name = {}
    member_gene_to_strain_name_dict = {}
    group_name_to_member_genes = {}
    next_group_id = 0

    for file in os.listdir(normalized_hits_dir):
        if not file.endswith('.m8'):
            continue

        logger.info(f'Processing file: {file}')
        hits_file_path = os.path.join(normalized_hits_dir, file)
        with open(hits_file_path) as f:
            for line in f:
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
        strain_to_members = {strain: ';'.join(sorted(members)) for strain, members in strain_to_members.items()}
        group_row_str = delimiter.join([f'OG_{group_index}'] + [strain_to_members.get(strain, '') for strain in sorted_strains])
        result += group_row_str + '\n'

    with open(putative_orthologs_path, 'w') as f:
        f.write(result)

    logger.info(f'Wrote putative orthogroups table to {putative_orthologs_path}')


def construct_table(logger, normalized_hits_dir, putative_orthologs_path, delimiter):
    if not os.path.exists(putative_orthologs_path):
        cluster_genes_to_connected_components(logger, normalized_hits_dir, putative_orthologs_path, delimiter)
    else:
        logger.info(f'Putative orthogroups table already exists at {putative_orthologs_path}')

    orthogroups_df = pd.read_csv(putative_orthologs_path, delimiter=delimiter)
    if len(orthogroups_df) == 0:
        raise ValueError('No putative ortholog groups were detected in your dataset. Please try to lower the '
                         'similarity parameters (see Advanced Options in the submission page) and re-submit your job.')

    orthogroups_df = orthogroups_df.sort_values(by=list(orthogroups_df.columns[1:])).reset_index(drop=True)
    orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(orthogroups_df.index))]
    orthogroups_df.to_csv(putative_orthologs_path, index=False)
    logger.info(f'Wrote sorted putative orthogroups table to {putative_orthologs_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('normalized_hits_dir',
                        help='path to a dir with all the hits files')
    parser.add_argument('putative_orthologs_path',
                        help='path to an output file in which the putative orthologs table will be written')
    parser.add_argument('--delimiter', help='delimiter for the input and output files', default=consts.CSV_DELIMITER)
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        construct_table(logger, args.normalized_hits_dir, args.putative_orthologs_path, args.delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
