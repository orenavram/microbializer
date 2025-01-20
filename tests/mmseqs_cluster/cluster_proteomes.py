import os
import subprocess
import sys
from sys import argv
import argparse
import traceback
import pandas as pd
import shutil
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger

MMSEQS_CLUSTER_MIN_SEQ_ID = 5
MMSEQS_CLUSTER_MIN_COVERAGE = 10

MMSEQS_CLUSTER_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1
MMSEQS_CLUSTER_REQUIRED_MEMORY_GB = '64'


def prepare_proteomes_subsets(logger, all_proteins_fasta_path, output_dir, clusters_file_path, min_seq_identity,
                              min_coverage, threads, num_of_clusters_in_orthogroup_inference):
    temp_outputs = os.path.join(output_dir, 'temp')
    os.makedirs(temp_outputs, exist_ok=True)

    # cluster all proteins
    cluster_result_prefix = os.path.join(temp_outputs, 'clusterRes')
    cmd = f'mmseqs easy-cluster {all_proteins_fasta_path} {cluster_result_prefix} {temp_outputs} --min-seq-id {min_seq_identity / 100} -c {min_coverage / 100} --cov-mode 0 --threads {threads} --remove-tmp-files'
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)

    # number clusters
    clusters_df = pd.read_csv(f'{cluster_result_prefix}_cluster.tsv', sep='\t',
                              names=['cluster_representative', 'cluster_member'])
    clusters_df['rep_id'] = pd.factorize(clusters_df['cluster_representative'])[0]

    # Create a mapping from original IDs to the target number of IDs using pandas.cut
    # Create 10 equally spaced bins and assign each ID to one of these bins
    clusters_df['cluster_id'] = pd.cut(clusters_df['rep_id'], bins=num_of_clusters_in_orthogroup_inference,
                                       labels=range(0, num_of_clusters_in_orthogroup_inference))

    # Convert to integer type
    clusters_df['cluster_id'] = clusters_df['cluster_id'].astype(int)

    clusters_df.to_csv(clusters_file_path, index=False)

    # remove temp dir
    shutil.rmtree(temp_outputs)

    # count number of sequence comparisons to be done
    cluster_counts = clusters_df.groupby('cluster_id').size().reset_index(name='count')
    number_of_comparison = sum(cluster_counts['count'] * (cluster_counts['count'] - 1) / 2)
    logger.info(f'Number of clusters: {len(cluster_counts)}, Number of sequence comparisons to be done: {number_of_comparison}')

    # write cluster members to separate files
    for cluster_id, group in clusters_df.groupby('cluster_id'):
        filename = os.path.join(output_dir, f"cluster_{cluster_id}.txt")
        group['cluster_member'].to_csv(filename, index=False, header=False)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('all_proteins_fasta_path', help='')
    parser.add_argument('output_dir', help='')
    parser.add_argument('clusters_file_path', help='')
    parser.add_argument('min_seq_identity', help='', type=float)
    parser.add_argument('min_coverage', help='', type=float)
    parser.add_argument('threads', help='')
    parser.add_argument('--num_of_clusters_in_orthogroup_inference', help='', default=5, type=int)
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir)

    logger.info(script_run_message)
    try:
        prepare_proteomes_subsets(logger, args.all_proteins_fasta_path, args.output_dir, args.clusters_file_path,
                                  args.min_seq_identity, args.min_coverage, args.threads,
                                  args.num_of_clusters_in_orthogroup_inference)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)