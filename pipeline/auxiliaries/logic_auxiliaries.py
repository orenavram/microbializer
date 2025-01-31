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
from Bio import SeqIO


from . import consts


def aggregate_mmseqs_scores(orthologs_scores_statistics_dir, paralogs_scores_statistics_dir, output_file, skip_paralogs):
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

    if not skip_paralogs:
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
    interval_length = math.ceil((end - start + 1) / number_of_intervals)

    # Create the intervals
    intervals = []
    i = start
    while i < end:
        interval_end = i + interval_length - 1
        intervals.append((i, interval_end))
        i = interval_end + 1

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


def convert_seq_identity_to_sensitivity(seq_identity):
    if seq_identity <= 0.3:
        sensitivity = 6
    elif seq_identity >= 0.8:
        sensitivity = 1.0
    else:
        sensitivity = 1.0 + (1.0 * (0.8 - seq_identity) * 10)

    return sensitivity


def fna_to_faa(logger, nucleotide_path, protein_path):
    with open(nucleotide_path, "r") as in_handle, open(protein_path, "w") as out_handle:
        # Iterate through each sequence record in the input file
        for record in SeqIO.parse(in_handle, "fasta"):
            # Translate the DNA sequence into a protein sequence
            translated_record = record.translate(id=True, name=True, description=True)

            # Write the translated record to the output file
            SeqIO.write(translated_record, out_handle, "fasta")

    logger.info(f'Translated fatsa file {nucleotide_path}. Output was written successfully to: {protein_path}')


def combine_orphan_genes_stats(orphan_genes_dir, output_dir):
    all_stat_dfs = []
    for file_name in os.listdir(orphan_genes_dir):
        if 'orphans_stats.csv' not in file_name:
            continue
        df = pd.read_csv(os.path.join(orphan_genes_dir, file_name), index_col=0)
        all_stat_dfs.append(df)

    combined_df = pd.concat(all_stat_dfs)
    combined_df.to_csv(os.path.join(output_dir, 'orphans_genes_stats.csv'))

    number_of_orphans_per_file = combined_df['Total orphans count'].to_dict()
    plot_genomes_histogram(number_of_orphans_per_file, output_dir, 'orphan_genes_count', 'Orphan genes count',
                           'Orphan genes count per Genome')


def split_ogs_to_jobs_inputs_files_by_og_sizes(orthogroups_df, step_tmp_dir, max_parallel_jobs):
    orthogroups_df['strains_count'] = orthogroups_df.notna().sum(axis=1) - 1
    orthogroups_df['paralogs_count'] = orthogroups_df.apply(lambda row:
                                                            sum(genes.count(';') for genes in row[1:-1] if pd.notna(genes)),
                                                            axis=1)
    orthogroups_df['genes_count'] = orthogroups_df['strains_count'] + orthogroups_df['paralogs_count']
    ogs_genes_count_df = orthogroups_df[['OG_name', 'genes_count']]
    ogs_genes_count_df.sort_values(by='genes_count', ascending=False, inplace=True)

    job_index_to_ogs = {i: [] for i in range(max_parallel_jobs)}
    job_index_to_genes_count = {i: 0 for i in range(max_parallel_jobs)}

    for _, row in ogs_genes_count_df.iterrows():
        # Find the job with the smallest current genes count
        job_index_with_min_genes_count = min(job_index_to_genes_count, key=job_index_to_genes_count.get)
        # Assign the current OG to this job
        job_index_to_ogs[job_index_with_min_genes_count].append(row["OG_name"])
        # Update the job's genes count
        job_index_to_genes_count[job_index_with_min_genes_count] += row["genes_count"]

    job_inputs_dir = os.path.join(step_tmp_dir, 'jobs_inputs')
    os.makedirs(job_inputs_dir, exist_ok=True)
    job_paths = []
    for job_index, ogs in job_index_to_ogs.items():
        job_path = os.path.join(job_inputs_dir, f'{job_index}.txt')
        with open(job_path, 'w') as f:
            f.write('\n'.join(map(str, ogs)))
        job_paths.append(job_path)

    # Remove the columns that were added
    orthogroups_df.drop(columns=['strains_count', 'paralogs_count', 'genes_count'], inplace=True)

    return job_paths
