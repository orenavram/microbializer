import sys
from sys import argv
import argparse
import logging
import os
import re
import time
import json

import matplotlib.pyplot as plt
from matplotlib import colors, patches as mpatches
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import get_job_logger, prepare_directories, submit_batch, wait_for_results
from calculate_cai import get_genome_to_W_vector

OG_NAME_PATTERN = re.compile(r'(OG_\d+)_cai.json')


def visualize_Ws_with_PCA(W_vectors, output_dir, logger):
    """
    output: 2D vectors (reduced with PCA) clustered with Kmeans
    """
    start_time = time.time()

    Ws_df = pd.DataFrame(W_vectors).transpose()
    Ws_values = Ws_df.values

    # Perform K-means clustering
    genome_count = len(W_vectors)
    if genome_count > 75:
        n_clusters = 5
    elif genome_count > 30:
        n_clusters = 4
    else:
        n_clusters = 3
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(Ws_values)

    # Perform PCA dimensionality reduction
    pca = PCA(n_components=2)
    Ws_values_reduced = pca.fit_transform(Ws_values)
    
    # Normalize cluster values for consistent colormap
    norm = colors.Normalize(vmin=0, vmax=n_clusters-1)

    # Create the scatter plot with color-coded clusters
    x = Ws_values_reduced[:, 0]
    y = Ws_values_reduced[:, 1]
    plt.axis('equal')
    scatter = plt.scatter(x, y, c=cluster_labels, cmap='viridis', alpha=0.4)
    plt.title("Relative Adaptiveness (W vectors)")
    
    cluster_legend = []
    for cluster in range(n_clusters):
        cluster_color = scatter.cmap(norm(cluster))
        legend_patch = mpatches.Patch(color=cluster_color, label=f'Cluster {cluster+1}')
        cluster_legend.append(legend_patch)
        
    plt.legend(handles=cluster_legend)

    # Save plot to output_dir
    plt.savefig(os.path.join(output_dir, 'Relative_Adaptiveness_scatter_plot.png'))
    plt.close()

    # Save point labels and coordinates to file
    point_labels_df = pd.DataFrame({'Genome': list(W_vectors.keys()), 'X': x, 'Y': y, 'Cluster': cluster_labels + 1})
    point_labels_df.to_csv(os.path.join(output_dir, 'point_labels.csv'), index=False)

    logger.info("Time for PCA:", time.time() - start_time)


def get_CAI_Data(cai_dir, output_dir, cai_table_path):
    """
    input:
    location of output dir to access individually calculated CAI data
   
    output: Histogram of CAI distribution and table of mean and standard
    deviation for each orthologous group
    """

    # Iterate through the directory of OG group CAIs and extract data
    CAI_Data = {}
    for filename in os.listdir(cai_dir):
        with open(os.path.join(cai_dir, filename), 'r') as cai_file:
            cai_info = json.load(cai_file)
            og_name = OG_NAME_PATTERN.match(filename).group(1)
            CAI_Data[og_name] = {"CAI_mean": cai_info["mean"], "CAI_std": cai_info["std"]}

    cai_df = pd.DataFrame(data=CAI_Data).transpose()
    cai_df.sort_values("mean", ascending=False, inplace=True)
    cai_df.to_csv(cai_table_path, index_label='OG_name')

    # Create histogram of CAI values
    plt.title("CAI distribution across OGs")
    plt.xlabel('CAI value')
    plt.ylabel('Frequency')
    plt.axis('auto')
    plt.hist(cai_df["mean"], bins=30)
    plt.savefig(os.path.join(output_dir, 'CAI_Histogram.png'))
    plt.close()


def analyze_codon_bias(ORF_dir, OG_dir, output_dir, cai_table_path, tmp_dir, src_dir, queue_name, error_file_path, logger):
    # 1. Calculate Ws
    step_number = '12_5_a'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_get_W'
    script_path = os.path.join(src_dir, 'steps/codon_bias/get_W.py')
    W_output_dir, W_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)

    all_cmds_params = []
    for orf_file_name in os.listdir(ORF_dir):
        if '02_ORFs' not in orf_file_name:  # Ignore system files that are automatically created sometimes
            continue
        orf_file_path = os.path.join(ORF_dir, orf_file_name)
        single_cmd_params = [orf_file_path,
                             W_output_dir,
                             W_tmp_dir]
        all_cmds_params.append(single_cmd_params)

    num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, W_tmp_dir,
                                               num_of_cmds_per_job=10,
                                               job_name_suffix='calc_W',
                                               queue_name=queue_name)

    # Passing logger also as times_logger since there is no convenient way here to get times_logger file path
    wait_for_results(logger, logger, os.path.split(script_path)[-1], W_tmp_dir,
                     num_of_batches, error_file_path)

    W_vectors = get_genome_to_W_vector(W_output_dir)
    if any(None in w_vector.values() for w_vector in W_vectors.values()):
        logger.info("At least one of the genomes has an incomplete W vector (might be because there were not "
                    "enough genes that were identified as HEGs). Codon bias analysis is therefore not possible.")
        return

    # 2. Make Graph
    visualize_Ws_with_PCA(W_vectors, output_dir, logger)

    # 3. Calculate CAIs
    step_number = '12_5_b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_calc_CAI'
    script_path = os.path.join(src_dir, 'steps/codon_bias/calculate_cai.py')
    cai_output_dir, cai_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)

    all_cmds_params = []
    number_of_ogs_per_job = 100
    start_og_index = 0
    max_og_index = len(os.listdir(args.OG_dir)) - 1
    while start_og_index <= max_og_index:
        stop_og_index = min(start_og_index + number_of_ogs_per_job - 1, max_og_index)
        single_cmd_params = [OG_dir,
                             W_output_dir,
                             start_og_index,
                             stop_og_index,
                             cai_output_dir]
        all_cmds_params.append(single_cmd_params)
        start_og_index += number_of_ogs_per_job

    num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, cai_tmp_dir,
                                               num_of_cmds_per_job=1,
                                               job_name_suffix='calc_cai',
                                               queue_name=queue_name)

    # Passing logger also as times_logger since there is no convenient way here to get times_logger file path
    wait_for_results(logger, logger, os.path.split(script_path)[-1], cai_tmp_dir,
                     num_of_batches, error_file_path)

    # 4. Make CAI table and histogram
    get_CAI_Data(cai_output_dir, output_dir, cai_table_path)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('ORF_dir', help='path to input fasta directory')
    parser.add_argument('OG_dir', help='path to input Orthologous group directory')
    parser.add_argument('output_dir', help='path to output directory')
    parser.add_argument('cai_table_path', help='path to the output CAI table of all OGs')
    parser.add_argument('tmp_dir', help='path to tmp directory')
    parser.add_argument('src_dir', help='path to pipeline directory')
    parser.add_argument('queue_name', help='queue to submit jobs to')
    parser.add_argument('error_file_path', help='path to error file')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        analyze_codon_bias(args.ORF_dir, args.OG_dir, args.output_dir, args.cai_table_path, args.tmp_dir, args.src_dir,
                           args.queue_name, args.error_file_path, logger)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')