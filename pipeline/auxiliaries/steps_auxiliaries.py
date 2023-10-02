import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def aggregate_ani_results(ani_output_dir):
    all_dfs = []
    for file_name in os.listdir(ani_output_dir):
        if not file_name.endswith('.tsv'):
            continue
        file_path = os.path.join(ani_output_dir, file_name)
        df = pd.read_csv(file_path, delimiter='\t',
                         names=['query', 'subject', 'ani_value', 'orthologous_segments', 'total_segments'])
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)

    query_short = []
    subject_short = []
    for index, row in combined_df.iterrows():
        query_short.append(os.path.splitext(os.path.basename(row['query']))[0])
        subject_short.append(os.path.splitext(os.path.basename(row['subject']))[0])

    combined_df['query'] = query_short
    combined_df['subject'] = subject_short

    ani_values_df = combined_df.pivot_table(index='query', columns='subject', values='ani_value')

    num_of_genomes = len(all_dfs)
    sns.set_context('paper', font_scale=1.4)
    if num_of_genomes <= 10:
        plt.subplots(figsize=(12, 10))
        plot = sns.heatmap(ani_values_df, annot=True, fmt='.1f', cmap='Blues')
    elif num_of_genomes <= 20:
        plt.subplots(figsize=(26, 20))
        plot = sns.heatmap(ani_values_df, cmap='Blues')
    else:
        plt.subplots(figsize=(53, 40))
        plot = sns.heatmap(ani_values_df, cmap='Blues')

    plt.tight_layout()
    fig = plot.get_figure()
    fig.savefig(os.path.join(ani_output_dir, 'heatmap.png'))

    ani_values_without_diagonal_df = ani_values_df.replace({100: np.nan})
    max_values = ani_values_without_diagonal_df.max(axis=1)
    max_columns = ani_values_without_diagonal_df.idxmax(axis=1)
    ani_values_df = ani_values_df.assign(max_value=max_values.values)
    ani_values_df = ani_values_df.assign(max_column=max_columns.values)
    ani_values_df.to_csv(os.path.join(ani_output_dir, 'ani_pairwise_values.csv'))


def mimic_prodigal_output(orfs_dir, output_orf_file_extension):
    for file_name in os.listdir(orfs_dir):
        file_path = os.path.join(orfs_dir, file_name)

        # edit headers of ORFs to match the structure of prodigal output
        fixed_content = ''
        with open(file_path, 'r') as orfs_file:
            for line in orfs_file:
                if line.startswith('>'):
                    fixed_content += f'{line.strip()} # START # END # 1 # \n'
                else:
                    fixed_content += line

        # override the old file with the fixed content
        with open(file_path, 'w') as f:
            f.write(fixed_content)

        # change file name to match the output of step 2
        os.rename(file_path, f'{os.path.splitext(file_path)[0]}.{output_orf_file_extension}')

