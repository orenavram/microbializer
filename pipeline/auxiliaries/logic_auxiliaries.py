import os
import re
import json

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import numpy as np

from . import consts


def aggregate_ani_results(ani_tmp_files, ani_output_dir):
    all_dfs = []
    for file_name in os.listdir(ani_tmp_files):
        if not file_name.endswith('.tsv'):
            continue
        file_path = os.path.join(ani_tmp_files, file_name)
        df = pd.read_csv(file_path, delimiter='\t',
                         names=['query', 'subject', 'ani_value', 'orthologous_segments', 'total_segments'])
        df['query'] = df['query'].apply(lambda path: os.path.splitext(os.path.basename(path))[0])
        df['subject'] = df['subject'].apply(lambda path: os.path.splitext(os.path.basename(path))[0])
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    ani_values_df = combined_df.pivot_table(index='query', columns='subject', values='ani_value')

    num_of_genomes = len(all_dfs)
    sns.set_context('paper', font_scale=1.4)
    if num_of_genomes <= 10:
        plt.subplots(figsize=(12, 10))
        sns.clustermap(ani_values_df, annot=True, fmt='.1f', cmap='Blues')
    elif num_of_genomes <= 20:
        plt.subplots(figsize=(26, 20))
        sns.clustermap(ani_values_df, cmap='Blues')
    else:
        plt.subplots(figsize=(53, 40))
        sns.clustermap(ani_values_df, cmap='Blues')

    plt.tight_layout()
    plt.savefig(os.path.join(ani_output_dir, 'ani_clustermap.png'))
    plt.clf()

    # Iterate over rows and find max value ignoring diagonal
    max_values = []
    max_values_columns = []
    for i, row in ani_values_df.iterrows():
        row_values = row.drop(i)  # Drop diagonal element
        max_values.append(row_values.max())
        max_values_columns.append(row_values.idxmax())

    ani_values_df['max_value'] = max_values
    ani_values_df['max_column'] = max_values_columns

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


def aggregate_mmseqs_scores(scores_statistics_dir, output_file):
    scores_means_per_strains_pair = {}
    scores_total_sum = 0
    scores_total_records = 0
    for scores_statistics_file in os.listdir(scores_statistics_dir):
        strains_names = os.path.splitext(scores_statistics_file)[0]
        with open(os.path.join(scores_statistics_dir, scores_statistics_file)) as fp:
            strains_statistics = json.load(fp)
        scores_means_per_strains_pair[strains_names] = strains_statistics['mean']
        scores_total_sum += strains_statistics['sum']
        scores_total_records += strains_statistics['number of records']

    scores_total_mean = scores_total_sum / scores_total_records

    scores_normalize_coefficients = {strains_names: strains_scores_mean / scores_total_mean
                                     for strains_names, strains_scores_mean in scores_means_per_strains_pair.items()}
    scores_statistics = {'mean_per_strain_pair': scores_means_per_strains_pair, 'total_scores_mean': scores_total_mean,
                         'scores_normalize_coefficients': scores_normalize_coefficients}
    with open(output_file, 'w') as fp:
        json.dump(scores_statistics, fp)

    return scores_normalize_coefficients


def get_strain_name(gene_name):
    return gene_name.split(':')[0]


def convert_required_sequence_identity_to_mmseqs_threshold(required_sequence_identity):
    # Formula taken from: https://github.com/soedinglab/MMseqs2/issues/777

    if required_sequence_identity <= 0.3:
        sens = 6
    elif required_sequence_identity > 0.8:
        sens = 1.0
    else:
        sens = 1.0 + (1.0 * (0.8 - required_sequence_identity) * 10)

    return min(sens + 1, 6)


def max_with_nan(x, y):
    if np.isnan(x):
        return y
    if np.isnan(y):
        return x
    return max(x, y)


def add_score_column_to_mmseqs_output(mmseqs_output_df):
    if consts.SIMILARITY_SCORE_CRITERION == consts.SimilarityScore.BITS:
        mmseqs_output_df['score'] = mmseqs_output_df['bits']
    else:  # consts.SIMILARITY_SCORE_CRITERION == consts.SimilarityScore.EVALUE
        mmseqs_output_df['score'] = -np.log10(mmseqs_output_df['evalue'])

        # Change Infinity scores (evalue = 0) to the max hit score + 1
        max_score = max(set(mmseqs_output_df['score']) - {np.inf})
        mmseqs_output_df.loc[mmseqs_output_df['score'] == np.inf, 'score'] = max_score + 1


def remove_bootstrap_values(in_tree_path, out_tree_path):
    with open(in_tree_path) as f:
        tree_as_str = f.read()

    tree_as_str = re.sub('\)\d+:', '):', tree_as_str)
    with open(out_tree_path, 'w') as f:
        f.write(tree_as_str)


def flatten(l):
    return [item.strip() for sublist in l for item in sublist if pd.notna(item)]


def get_all_genes_in_table(df):
    # Expect as input a dataframe with genes in the cells (without an OG column)
    all_df_values = flatten(df.values)

    all_genes = flatten([value.split(';') for value in all_df_values])
    return all_genes


def plot_genomes_histogram(data, output_dir, output_file_name, title, xlabel):
    # data is expected to be: {'genome1': 54, 'genome2': 20, ...}

    with open(os.path.join(output_dir, f'{output_file_name}.json'), 'w') as fp:
        json.dump(data, fp)

    output_df = pd.DataFrame.from_dict(data, orient='index', columns=[title])
    output_df.index.name = 'Genome'
    output_df.to_csv(os.path.join(output_dir, f'{output_file_name}.csv'))

    sns.histplot(output_df, x=title, kde=True)
    plt.title(f'Distribution of {title}')
    plt.xlabel(xlabel)
    plt.ylabel('Genomes count')
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))  # make y-axis integer
    plt.savefig(os.path.join(output_dir, f'{output_file_name}.png'))

    plt.clf()
