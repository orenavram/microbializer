import os
import subprocess
import sys
from sys import argv
import argparse
import logging
import pandas as pd
import shutil
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def prepare_proteomes_subsets(logger, translated_orfs_dir, output_dir, clusters_file_path, min_seq_identity,
                              min_coverage, threads, do_cluster):
    temp_outputs = os.path.join(output_dir, 'temp')
    os.makedirs(temp_outputs, exist_ok=True)

    # create a fasta files of all proteomes
    all_proteins_fasta_path = os.path.join(output_dir, 'all_proteomes.faa')
    cmd = f"cat {os.path.join(translated_orfs_dir, '*')} > {all_proteins_fasta_path}"
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)

    if do_cluster:
        # cluster all proteins
        cluster_result_prefix = os.path.join(temp_outputs, 'clusterRes')
        cmd = f'mmseqs easy-cluster {all_proteins_fasta_path} {cluster_result_prefix} {temp_outputs} --min-seq-id {min_seq_identity / 100} -c {min_coverage / 100} --cov-mode 0 --threads {threads} --remove-tmp-files'
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)

        # number clusters
        clusters_df = pd.read_csv(f'{cluster_result_prefix}_cluster.tsv', sep='\t',
                                  names=['cluster_representative', 'cluster_member'])
        clusters_df['rep_id'] = pd.factorize(clusters_df['cluster_representative'])[0]

        # Number of unique IDs to reduce to
        target_num_ids = 5

        # Create a mapping from original IDs to the target number of IDs using pandas.cut
        # Create 10 equally spaced bins and assign each ID to one of these bins
        clusters_df['cluster_id'] = pd.cut(clusters_df['rep_id'], bins=target_num_ids, labels=range(0, target_num_ids))

        # Convert to integer type
        clusters_df['cluster_id'] = clusters_df['cluster_id'].astype(int)
    else:
        all_record_ids = [record.id for record in SeqIO.parse(all_proteins_fasta_path, "fasta")]
        clusters_df = pd.DataFrame({
            'cluster_id': 0,
            'cluster_representative': 'nan',
            'cluster_member': all_record_ids
        })

    clusters_df.to_csv(clusters_file_path, index=False)

    # remove temp dir
    shutil.rmtree(temp_outputs)

    # count number of sequence comparisons to be done
    cluster_counts = clusters_df.groupby('cluster_id').size().reset_index(name='count')
    number_of_comparison = sum(cluster_counts['count'] * (cluster_counts['count'] - 1) / 2)
    logger.info(f'Number of clusters: {len(cluster_counts)}, Number of sequence comparisons to be done: {number_of_comparison}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('translated_orfs_dir', help='')
    parser.add_argument('output_dir', help='')
    parser.add_argument('clusters_file_path', help='')
    parser.add_argument('min_seq_identity', help='', type=float)
    parser.add_argument('min_coverage', help='', type=float)
    parser.add_argument('threads', help='')
    parser.add_argument('--do_cluster', help='', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir)

    logger.info(script_run_message)
    try:
        prepare_proteomes_subsets(logger, args.translated_orfs_dir, args.output_dir, args.clusters_file_path,
                                  args.min_seq_identity, args.min_coverage, args.threads, args.do_cluster)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            f.write(f'Internal Error in {__file__}: {e}\n')
