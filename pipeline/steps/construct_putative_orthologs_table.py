import argparse
from pathlib import Path
import sys
from collections import defaultdict
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step


def cluster_genes_to_connected_components(logger, normalized_hits_dir, genome_names_file_path, putative_orthologs_path):
    member_gene_to_group_name = {}
    member_gene_to_strain_name_dict = {}
    group_name_to_member_genes = {}
    next_group_id = 0

    for hits_file_path in normalized_hits_dir.glob('*.m8'):
        logger.info(f'Processing file: {hits_file_path}')
        with open(hits_file_path) as f:
            # Read header
            line = f.readline()
            strain1, strain2 = line.rstrip().split(',')[:2]

            # Read the rest of the file
            for line in f:
                line_tokens = line.rstrip().split(',')
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
    with open(genome_names_file_path) as genome_names_fp:
        strain_names = [line.strip() for line in genome_names_fp if line.strip()]
    sorted_strains = sorted(strain_names)

    header = ','.join(['OG_name'] + sorted_strains) + '\n'
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
        group_row_str = ','.join(
            [f'OG_{group_index}'] + [strain_to_members.get(strain, '') for strain in sorted_strains])
        result += group_row_str + '\n'

    with open(putative_orthologs_path, 'w') as f:
        f.write(result)

    logger.info(f'Wrote putative orthogroups table to {putative_orthologs_path}')


def construct_table(logger, normalized_hits_dir, genome_names_file_path, putative_orthologs_path):
    if not putative_orthologs_path.exists():
        cluster_genes_to_connected_components(logger, normalized_hits_dir, genome_names_file_path, putative_orthologs_path)
    else:
        logger.info(f'Putative orthogroups table already exists at {putative_orthologs_path}')

    orthogroups_df = pd.read_csv(putative_orthologs_path, dtype=str, engine='pyarrow', dtype_backend='pyarrow')
    if len(orthogroups_df) == 0:
        raise ValueError('No putative ortholog groups were detected in your dataset. Please try to lower the '
                         'similarity parameters (see Advanced Options in the submission page) and re-submit your job.')

    orthogroups_df = orthogroups_df.sort_values(by=list(orthogroups_df.columns[1:])).reset_index(drop=True)
    orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(orthogroups_df.index))]
    orthogroups_df.to_csv(putative_orthologs_path, index=False)
    logger.info(f'Wrote sorted putative orthogroups table to {putative_orthologs_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('normalized_hits_dir', type=Path,
                        help='path to a dir with all the hits files')
    parser.add_argument('genome_names_file_path', type=Path)
    parser.add_argument('putative_orthologs_path', type=Path,
                        help='path to an output file in which the putative orthologs table will be written')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, construct_table, args.normalized_hits_dir, args.genome_names_file_path, args.putative_orthologs_path)
