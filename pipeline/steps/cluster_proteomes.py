import os
import subprocess
import sys
from sys import argv
import argparse
import logging
import pandas as pd
import shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def prepare_proteomes_subsets(logger, translated_orfs_dir, output_dir, clusters_file_path, min_seq_identity, min_coverage, threads):
    temp_outputs = os.path.join(output_dir, 'temp')
    os.makedirs(temp_outputs, exist_ok=True)

    # create a fasta files of all proteomes
    all_proteins_fasta_path = os.path.join(temp_outputs, 'all_proteomes.faa')
    cmd = f"cat {os.path.join(translated_orfs_dir, '*')} > {all_proteins_fasta_path}"
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)

    # cluster all proteins
    cluster_result_prefix = os.path.join(temp_outputs, 'clusterRes')
    cluster_tmp_dir = os.path.join(output_dir, 'tmp')
    cmd = f'mmseqs easy-cluster {all_proteins_fasta_path} {cluster_result_prefix} {cluster_tmp_dir} --min-seq-id {min_seq_identity} -c {min_coverage} --cov-mode 0 --threads {threads} --remove-tmp-files'
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)

    # number clusters
    clusters_df = pd.read_csv(f'{cluster_result_prefix}_cluster.tsv', sep='\t',
                              names=['cluster_representative', 'cluster_member'])
    clusters_df['cluster_id'] = pd.factorize(clusters_df['cluster_representative'])[0]
    clusters_df.to_csv(clusters_file_path, index=False)

    # remove temp dir
    shutil.rmtree(temp_outputs)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('translated_orfs_dir', help='')
    parser.add_argument('output_dir', help='')
    parser.add_argument('clusters_file_path', help='')
    parser.add_argument('min_seq_identity', help='')
    parser.add_argument('min_coverage', help='')
    parser.add_argument('threads', help='')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        prepare_proteomes_subsets(logger, args.translated_orfs_dir, args.output_dir, args.clusters_file_path,
                                  args.min_seq_identity, args.min_coverage, args.threads)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
