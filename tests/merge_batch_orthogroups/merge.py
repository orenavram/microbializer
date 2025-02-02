import os

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'output')

pseudo_orthogroups_file_path = os.path.join(SCRIPT_DIR, 'pseudo_orthogroups.csv')
sub_orthogroups_dir_path = os.path.join(SCRIPT_DIR, 'sub_orthogroups')


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pseudo_orthogroups_df = pd.read_csv(pseudo_orthogroups_file_path)

    for sub_orthogroups_file_name in os.listdir(sub_orthogroups_dir_path):
        batch_id = os.path.splitext(sub_orthogroups_file_name)[0].split('_')[-1]

        sub_orthogroups_file_path = os.path.join(sub_orthogroups_dir_path, sub_orthogroups_file_name)
        sub_orthogroups_df = pd.read_csv(sub_orthogroups_file_path)
        sub_orthogroups_df.drop(columns=['OG_name'], inplace=True)
        sub_orthogroups_df.set_index('representative_gene', inplace=True)

        # Fill rows of pseudo_orthogroups_df where pseudo_genome_{batch_id} contains only 1 gene
        pseudo_orthogroups_df = pseudo_orthogroups_df.merge(
            sub_orthogroups_df, how='left', left_on=f'pseudo_genome_{batch_id}', right_on='representative_gene')

        # Fill rows of pseudo_orthogroups_df where pseudo_genome_{batch_id} contains multiple genes
        pseudo_orthogroups_with_paralogs_in_pseudo_genome = pseudo_orthogroups_df[
            pseudo_orthogroups_df[f'pseudo_genome_{batch_id}'].str.contains(';', na=False)]
        for i, row in pseudo_orthogroups_with_paralogs_in_pseudo_genome.iterrows():
            pseudo_genes = row[f'pseudo_genome_{batch_id}'].split(';')
            row_sub_orthogroups = pd.concat([sub_orthogroups_df.loc[pseudo_gene] for pseudo_gene in pseudo_genes], axis=1)
            strain_to_genes = dict(
                row_sub_orthogroups.apply(lambda row: ';'.join(sorted(row.dropna().astype(str))), axis=1))
            for strain, genes in strain_to_genes.items():
                pseudo_orthogroups_df.at[i, strain] = genes if genes else np.nan

        pseudo_orthogroups_df.drop(columns=[f'pseudo_genome_{batch_id}'], inplace=True)

    pseudo_orthogroups_df = pseudo_orthogroups_df.sort_values(
        by=list(pseudo_orthogroups_df.columns[1:])).reset_index(drop=True)
    pseudo_orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(pseudo_orthogroups_df.index))]
    merged_orthogroups_file_path = os.path.join(OUTPUT_DIR, 'merged_orthogroups.csv')
    pseudo_orthogroups_df.to_csv(merged_orthogroups_file_path, index=False)


if __name__ == '__main__':
    main()
