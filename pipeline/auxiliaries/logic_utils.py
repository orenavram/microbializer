from __future__ import annotations

import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import numpy as np

from . import consts


def aggregate_mmseqs_scores(logger, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir, output_file,
                            skip_paralogs):
    scores_means_per_strains_pair = {}
    scores_total_sum = 0
    scores_total_records = 0
    for scores_statistics_file in orthologs_scores_statistics_dir.glob('*.stats'):
        strains_names = scores_statistics_file.stem
        with open(scores_statistics_file) as fp:
            strains_statistics = json.load(fp)
        scores_means_per_strains_pair[strains_names] = strains_statistics['mean']
        scores_total_sum += strains_statistics['sum']
        scores_total_records += strains_statistics['number of records']

    if not skip_paralogs:
        for scores_statistics_file in paralogs_scores_statistics_dir.glob('*.stats'):
            strains_names = scores_statistics_file.stem
            with open(scores_statistics_file) as fp:
                strains_statistics = json.load(fp)
            scores_means_per_strains_pair[strains_names] = strains_statistics['mean']
            scores_total_sum += strains_statistics['sum']
            scores_total_records += strains_statistics['number of records']

    if scores_total_records == 0:
        logger.warning('No orthologs or paralogs were found. No scores normalization will be done.')
        return None

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


def plot_genomes_histogram(data, output_dir, output_file_name, title, xlabel):
    # data is expected to be: {'genome1': 54, 'genome2': 20, ...}

    with open(output_dir / f'{output_file_name}.json', 'w') as fp:
        json.dump(data, fp)

    output_df = pd.DataFrame.from_dict(data, orient='index', columns=[title])
    output_df.index.name = 'Genome'
    output_df.to_csv(output_dir / f'{output_file_name}.csv')

    sns.histplot(output_df, x=title, kde=True)
    plt.title(f'Distribution of {title}', fontsize=25, loc='center', wrap=True)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel('Genomes count', fontsize=20)
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))  # make y-axis integer
    plt.tight_layout()
    plt.savefig(output_dir / f'{output_file_name}.png', dpi=600)

    plt.clf()


def convert_seq_identity_to_sensitivity(seq_identity):
    if seq_identity <= 0.3:
        sensitivity = 6
    elif seq_identity >= 0.8:
        sensitivity = 1.0
    else:
        sensitivity = 1.0 + (1.0 * (0.8 - seq_identity) * 10)

    return sensitivity


def combine_orphan_genes_stats(orphan_genes_dir, output_dir):
    all_stat_dfs = []
    for file_path in orphan_genes_dir.glob('*orphans_stats.csv'):
        df = pd.read_csv(file_path, index_col=0)
        all_stat_dfs.append(df)

    combined_df = pd.concat(all_stat_dfs)
    combined_df.to_csv(output_dir / 'orphans_genes_stats.csv')

    number_of_orphans_per_file = combined_df['Total orphans count'].to_dict()
    plot_genomes_histogram(number_of_orphans_per_file, output_dir, 'orphan_genes_count', 'Orphan genes count',
                           'Orphan genes count per Genome')


def count_strains_and_genes_in_ogs(orthogroups_df):
    orthogroups_df['strains_count'] = orthogroups_df.count(axis=1) - 1
    orthogroups_df['paralogs_count'] = orthogroups_df.apply(
        lambda row: sum(genes.count(';') for genes in row[1:-1] if pd.notna(genes)), axis=1)
    orthogroups_df['genes_count'] = orthogroups_df['strains_count'] + orthogroups_df['paralogs_count']

    return orthogroups_df[['OG_name', 'strains_count', 'genes_count']]


def split_ogs_to_jobs_inputs_files_by_og_sizes(orthogroups_df, step_tmp_dir, max_parallel_jobs):
    orthogroups_sizes_df = count_strains_and_genes_in_ogs(orthogroups_df)
    ogs_genes_count_df = orthogroups_sizes_df[['OG_name', 'genes_count']].sort_values(by='genes_count', ascending=False)

    job_index_to_ogs = {i: [] for i in range(max_parallel_jobs)}
    job_index_to_genes_count = {i: 0 for i in range(max_parallel_jobs)}

    for _, row in ogs_genes_count_df.iterrows():
        # Find the job with the smallest current genes count
        job_index_with_min_genes_count = min(job_index_to_genes_count, key=job_index_to_genes_count.get)
        # Assign the current OG to this job
        job_index_to_ogs[job_index_with_min_genes_count].append(row["OG_name"])
        # Update the job's genes count
        job_index_to_genes_count[job_index_with_min_genes_count] += row["genes_count"]

    job_inputs_dir = step_tmp_dir / 'jobs_inputs'
    job_inputs_dir.mkdir(parents=True, exist_ok=True)
    job_paths = []
    for job_index, ogs in job_index_to_ogs.items():
        job_path = job_inputs_dir / f'{job_index}.txt'
        with open(job_path, 'w') as f:
            f.write('\n'.join(map(str, ogs)))
        job_paths.append(job_path)

    # Remove the columns that were added
    orthogroups_df.drop(columns=['strains_count', 'paralogs_count', 'genes_count'], inplace=True)

    return job_paths


def count_and_plot_orthogroups_sizes(final_orthogroups_file_path, group_sizes_path):
    final_orthologs_table_df = pd.read_csv(final_orthogroups_file_path, dtype=str)
    orthogroups_sizes_df = count_strains_and_genes_in_ogs(final_orthologs_table_df)
    orthogroups_sizes_df = orthogroups_sizes_df.rename(columns={'strains_count': 'OG size (number of genomes)',
                                                                'genes_count': 'OG size (total number of genes)'})
    orthogroups_sizes_df.to_csv(group_sizes_path / 'groups_sizes.csv', index=False)

    group_sizes = orthogroups_sizes_df.set_index('OG_name')['OG size (number of genomes)']
    sns.histplot(x=group_sizes, discrete=True)
    if len(np.unique(group_sizes)) < 10:
        plt.xticks(np.unique(group_sizes))
    plt.title('Orthogroups sizes distribution', fontsize=25, loc='center', wrap=True)
    plt.xlabel('OG size (number of genomes)', fontsize=20)
    plt.ylabel('Count of OGs of each OG size', fontsize=20)
    plt.tight_layout()
    plt.savefig(group_sizes_path / 'groups_sizes.png', dpi=600)
    plt.clf()
