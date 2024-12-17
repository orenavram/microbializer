from __future__ import annotations

import os
import re
import json
import math

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import numpy as np


from . import consts


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


def aggregate_mmseqs_scores(orthologs_scores_statistics_dir, paralogs_scores_statistics_dir, output_file):
    scores_means_per_strains_pair = {}
    scores_total_sum = 0
    scores_total_records = 0
    for scores_statistics_file in os.listdir(orthologs_scores_statistics_dir):
        if not scores_statistics_file.endswith('.stats'):
            continue
        strains_names = os.path.splitext(scores_statistics_file)[0]
        with open(os.path.join(orthologs_scores_statistics_dir, scores_statistics_file)) as fp:
            strains_statistics = json.load(fp)
        scores_means_per_strains_pair[strains_names] = strains_statistics['mean']
        scores_total_sum += strains_statistics['sum']
        scores_total_records += strains_statistics['number of records']

    for scores_statistics_file in os.listdir(paralogs_scores_statistics_dir):
        if not scores_statistics_file.endswith('.stats'):
            continue
        strains_names = os.path.splitext(scores_statistics_file)[0]
        with open(os.path.join(paralogs_scores_statistics_dir, scores_statistics_file)) as fp:
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


def max_with_nan(x, y):
    if np.isnan(x):
        return y
    if np.isnan(y):
        return x
    return max(x, y)


def add_score_column_to_mmseqs_output(mmseqs_output_df):
    if consts.SIMILARITY_SCORE_CRITERION == consts.SimilarityScore.BITS:
        mmseqs_output_df['score'] = mmseqs_output_df['bits']
    elif consts.SIMILARITY_SCORE_CRITERION == consts.SimilarityScore.EVALUE:
        mmseqs_output_df['score'] = -np.log10(mmseqs_output_df['evalue'])

        # Change Infinity scores (evalue = 0) to the max hit score + 1
        max_score = max(set(mmseqs_output_df['score']) - {np.inf})
        mmseqs_output_df.loc[mmseqs_output_df['score'] == np.inf, 'score'] = max_score + 1
    else:
        raise ValueError(f"Unknown similarity score criterion: {consts.SIMILARITY_SCORE_CRITERION}")


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
    plt.title(f'Distribution of {title}', fontsize=20, loc='center', wrap=True)
    plt.xlabel(xlabel, fontsize=15)
    plt.ylabel('Genomes count', fontsize=15)
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))  # make y-axis integer
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{output_file_name}.png'), dpi=600)

    plt.clf()


def update_progressbar(progressbar_file_path, step_name_finished):
    df = pd.read_csv(progressbar_file_path)
    df.loc[df['Step'] == step_name_finished, 'Finished'] = True
    df.to_csv(progressbar_file_path, index=False)


def define_intervals(start, end, number_of_intervals):
    # Calculate the interval length
    interval_length = math.ceil((end - start) / number_of_intervals)

    # Create the intervals
    intervals = []
    i = start
    while i < end:
        interval_end = i + interval_length
        intervals.append((i, interval_end))
        i = interval_end

    # Adjust the last interval to ensure it ends exactly at end
    intervals[-1] = (intervals[-1][0], end)

    return intervals


def get_directory_size_in_gb(directory):
    total_size = 0
    for file in os.listdir(directory):
        file_path = os.path.join(directory, file)
        # Skip broken symlinks and non-files
        if os.path.isfile(file_path):
            total_size += os.path.getsize(file_path)
    # Convert bytes to gigabytes and round up
    size_in_gb = total_size / (1024 ** 3)
    return size_in_gb
