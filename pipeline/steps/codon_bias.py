import sys
from sys import argv
import argparse
import logging
import os
import time
import json
import subprocess
import shutil
from multiprocessing import Pool
import re
import traceback

import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import CodonAdaptationIndex
import matplotlib.pyplot as plt
from matplotlib import colors, patches as mpatches
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries.logic_auxiliaries import get_strain_name
from auxiliaries import consts

MAX_HITS_TO_KEEP_FOR_EACH_REFERENCE_HEG = 3
BLAST_IDENTITY_PERCENT_THRESHOLD = 35
BLAST_EVALUE_THRESHOLD = 0.01


def clean_seq_records(records):
    """
    Remove all non-complete codons from the records (e.g. NGT)
    """
    records_were_cleaned = False

    for record in records:
        clean_seq = ''
        for i in range(0, len(record.seq), 3):
            codon = record.seq[i: i + 3]
            if re.compile(r'^[ATGC][ATGC][ATGC]$', re.IGNORECASE).match(str(codon)):
                clean_seq += codon
        if record.seq != clean_seq:
            records_were_cleaned = True
        record.seq = clean_seq

    return records_were_cleaned


def find_HEGs_in_orf_file(ORFs_file, genome_name, hegs_output_dir, logger):
    # Create blast db from ORF file
    db_dir = os.path.join(hegs_output_dir, f'{genome_name}_db')
    db_name = os.path.join(db_dir, genome_name + '.db')
    cmd = f'makeblastdb -in {ORFs_file} -out {db_name} -dbtype nucl'
    logger.info('Making blastdb with command: \n' + cmd)
    subprocess.run(cmd, shell=True)

    # Query blast db with ecoli HEGs reference file
    hegs_hits_file = os.path.join(hegs_output_dir, genome_name + '_HEG_hits.tsv')
    cmd = f'tblastn -db {db_name} -query {consts.HEGS_ECOLI_FILE_PATH} -out {hegs_hits_file} -outfmt 6 ' \
          f'-max_target_seqs {MAX_HITS_TO_KEEP_FOR_EACH_REFERENCE_HEG}'
    logger.info("Finding Hits with command: \n" + cmd)
    subprocess.run(cmd, shell=True)

    shutil.rmtree(db_dir, ignore_errors=True)

    # Filter hits to find actual HEGs and write their names into a file
    hegs_df = pd.read_csv(hegs_hits_file, delimiter='\t', names=consts.BLAST_OUTPUT_HEADER)
    hegs_df_filtered = hegs_df.loc[(hegs_df['identity_percent'] > BLAST_IDENTITY_PERCENT_THRESHOLD) &
                                   (hegs_df['evalue'] < BLAST_EVALUE_THRESHOLD)]
    hegs_names = set(hegs_df_filtered['subject'])

    return hegs_names


def get_W(ORFs_file, hegs_output_dir, logger):
    logger.info("Getting HEGs from the ORFs file: " + str(ORFs_file))

    genome_name = os.path.splitext(os.path.basename(ORFs_file))[0]

    # Identify HEGs in the ORFs file
    HEGs_names = find_HEGs_in_orf_file(ORFs_file, genome_name, hegs_output_dir, logger)
    HEGs_records = [record for record in SeqIO.parse(ORFs_file, "fasta") if record.id in HEGs_names]
    SeqIO.write(HEGs_records, os.path.join(hegs_output_dir, genome_name + "_HEGs.fa"), "fasta")  # Write HEGs to a file for logging and debugging

    # Find genome codon index of HEGs
    records_were_cleaned = clean_seq_records(HEGs_records)
    if records_were_cleaned:
        logger.warning(f'Non-complete or illegal codons were found in the ORFs of {genome_name}. They were removed.')
        SeqIO.write(HEGs_records, os.path.join(hegs_output_dir, genome_name + "_HEGs_cleaned.fa"), "fasta")  # Write cleaned HEGs to a file for logging and debugging
    genome_codon_index = CodonAdaptationIndex(HEGs_records)

    return genome_name, genome_codon_index


def visualize_Ws_with_PCA(W_vectors, output_dir, logger):
    """
    output: 2D vectors (reduced with PCA) clustered with Kmeans
    """
    start_time = time.time()

    # Perform K-means clustering
    genome_count = len(W_vectors)
    if genome_count > 75:
        n_clusters = 5
    elif genome_count > 30:
        n_clusters = 4
    else:
        n_clusters = 3
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(W_vectors.values)

    # Perform PCA dimensionality reduction
    pca = PCA(n_components=2)
    Ws_values_reduced = pca.fit_transform(W_vectors.values)
    explained_variance_ratio = pca.explained_variance_ratio_ * 100

    # Normalize cluster values for consistent colormap
    norm = colors.Normalize(vmin=0, vmax=n_clusters - 1)

    # Create the scatter plot with color-coded clusters
    x = Ws_values_reduced[:, 0]
    y = Ws_values_reduced[:, 1]
    plt.axis('equal')
    scatter = plt.scatter(x, y, c=cluster_labels, cmap='viridis', alpha=0.4)
    plt.title("Relative Adaptiveness (W vectors) of Genomes", fontsize=20, loc='center', wrap=True)
    plt.xlabel(f"PC1 ({explained_variance_ratio[0]:.1f}%)", fontsize=15)
    plt.ylabel(f"PC2 ({explained_variance_ratio[1]:.1f}%)", fontsize=15)

    cluster_legend = []
    for cluster in range(n_clusters):
        cluster_color = scatter.cmap(norm(cluster))
        legend_patch = mpatches.Patch(color=cluster_color, label=f'Cluster {cluster + 1}')
        cluster_legend.append(legend_patch)

    plt.legend(handles=cluster_legend)

    # Save plot to output_dir
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Relative_Adaptiveness_scatter_plot.png'), dpi=600)
    plt.close()

    # Save point labels and coordinates to file
    point_labels_df = pd.DataFrame({'Genome': W_vectors.index, 'X': x, 'Y': y, 'Cluster': cluster_labels + 1})
    point_labels_df.to_csv(os.path.join(output_dir, 'Relative_Adaptiveness_scatter_plot_clusters.csv'), index=False)

    logger.info("Time for PCA:", time.time() - start_time)


def calculate_cai(OG_dir, OG_index, genomes_codon_indexes, cais_output_dir):
    OG_path = os.path.join(OG_dir, f'OG_{OG_index}_dna.fas')
    output_file_path = os.path.join(cais_output_dir, f'OG_{OG_index}.json')

    logger.info(f'Calculating CAI for genes in OG {OG_path}. Will save the results in {output_file_path}')

    cai_info = {}
    OG_records = list(SeqIO.parse(OG_path, 'fasta'))
    clean_seq_records(OG_records)
    for record in OG_records:
        genome_name = get_strain_name(record.id)
        cai_info[record.id] = genomes_codon_indexes[genome_name].calculate(record)

    cai_values = list(cai_info.values())
    cai_info['CAI_mean'] = np.mean(cai_values)
    cai_info['CAI_std'] = np.std(cai_values)

    with open(output_file_path, 'w') as output_file:
        json.dump(cai_info, output_file)

    summarized_cai_info = {key: value for key, value in cai_info.items() if key in ['CAI_mean', 'CAI_std']}
    return f'OG_{OG_index}', summarized_cai_info


def plot_CAI_histogram(ogs_cai_info_df, output_dir):
    plt.title("CAI Mean distribution across OGs", fontsize=20, loc='center', wrap=True)
    plt.xlabel('CAI value', fontsize=15)
    plt.ylabel('Frequency', fontsize=15)
    plt.axis('auto')
    plt.hist(ogs_cai_info_df["CAI_mean"], bins=30)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'CAI_histogram.png'), dpi=600)
    plt.close()


def analyze_codon_bias(ORF_dir, OG_dir, output_dir, cai_table_path, tmp_dir, cpus, logger):
    # 1. Calculate W vector for each genome
    hegs_output_dir = os.path.join(tmp_dir, 'HEGs')
    os.makedirs(hegs_output_dir, exist_ok=True)
    with Pool(processes=cpus) as pool:
        params = [(os.path.join(ORF_dir, orf_file), hegs_output_dir, logger) for orf_file in os.listdir(ORF_dir)]
        genomes_codon_indexes = pool.starmap(get_W, params)

    genomes_codon_indexes = dict(genomes_codon_indexes)
    genomes_W_vectors_df = pd.DataFrame.from_dict(genomes_codon_indexes, orient='index')
    genomes_W_vectors_df.to_csv(os.path.join(output_dir, 'W_vectors.csv'), index_label='Genome')

    # 2. Visualize W vectors with PCA
    if len(genomes_W_vectors_df) >= 4:
        visualize_Ws_with_PCA(genomes_W_vectors_df, output_dir, logger)

    # 3. Calculate CAI for each OG
    cais_output_dir = os.path.join(tmp_dir, 'OGs_CAIs')
    os.makedirs(cais_output_dir, exist_ok=True)
    with Pool(processes=cpus) as pool:
        params = [(OG_dir, OG_index, genomes_codon_indexes, cais_output_dir) for OG_index in range(len(os.listdir(OG_dir)))]
        ogs_cai_info = pool.starmap(calculate_cai, params)

    ogs_cai_info_df = pd.DataFrame.from_dict(dict(ogs_cai_info), orient='index')
    ogs_cai_info_df.sort_values("CAI_mean", ascending=False, inplace=True)
    ogs_cai_info_df.to_csv(cai_table_path, index_label='OG_name')

    # 4. Aggregate CAI data
    plot_CAI_histogram(ogs_cai_info_df, output_dir)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('ORF_dir', help='path to input fasta directory')
    parser.add_argument('OG_dir', help='path to input Orthologous group directory')
    parser.add_argument('output_dir', help='path to output directory')
    parser.add_argument('cai_table_path', help='path to the output CAI table of all OGs')
    parser.add_argument('tmp_dir', help='path to tmp directory')
    parser.add_argument('cpus', help='number of cpus to use')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        analyze_codon_bias(args.ORF_dir, args.OG_dir, args.output_dir, args.cai_table_path, args.tmp_dir, int(args.cpus), logger)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
