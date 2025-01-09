from sys import argv
import argparse
import logging
import os
import sys
import pandas as pd
from collections import defaultdict
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries.logic_auxiliaries import get_strain_name


def get_verified_clusters_set(verified_clusters_path):
    return set([os.path.splitext(file)[0]
                for file in os.listdir(verified_clusters_path) if file.endswith('verified_cluster')])


def finalize_table(logger, putative_orthologs_path, verified_clusters_path, finalized_table_path):
    putative_orthologs_df = pd.read_csv(putative_orthologs_path)

    # Keep only the verified clusters
    verified_clusters_set = get_verified_clusters_set(verified_clusters_path)
    verified_clusters_df = putative_orthologs_df.loc[putative_orthologs_df['OG_name'].isin(verified_clusters_set)]

    logger.info(f'Found {len(verified_clusters_df.index)} verified clusters from previous steps, out of '
                f'{len(putative_orthologs_df.index)} putative clusters overall. Starting to add putative clusters '
                f'that were split by MCL algorithm...')

    # Iterate through new clusters (that were split from the putative clusters) and create Series objects from them
    new_clusters = []
    for cluster_file_name in os.listdir(verified_clusters_path):
        if cluster_file_name.endswith('verified_cluster'):
            continue
        with open(os.path.join(verified_clusters_path, cluster_file_name), 'r') as new_cluster:
            genes = new_cluster.readline().strip().split('\t')
        strain_to_genes = defaultdict(list)
        for gene in genes:
            strain = get_strain_name(gene)
            strain_to_genes[strain].append(gene)
        strain_to_genes = {strain: ';'.join(sorted(genes)) for strain, genes in strain_to_genes.items()}
        strain_to_genes['OG_name'] = ''  # Set a temp empty OG name
        new_cluster = pd.Series(strain_to_genes)
        new_clusters.append(new_cluster)

    # Merge new clusters and verified clusters
    new_clusters_df = pd.DataFrame(new_clusters)
    all_clusters_df = pd.concat([verified_clusters_df, new_clusters_df], ignore_index=True)
    all_clusters_df = all_clusters_df.sort_values(by=list(all_clusters_df.columns[1:])).reset_index(drop=True)
    all_clusters_df['OG_name'] = [f'OG_{i}' for i in range(len(all_clusters_df.index))]
    all_clusters_df.to_csv(finalized_table_path, index=False)

    logger.info(f'Finished adding split clusters. OG table now contains {len(all_clusters_df.index)} groups, '
                f'and was written to {finalized_table_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('putative_orthologs_path', help='path to a file with the putative orthologs sets')
    parser.add_argument('verified_clusters_path', help='path to a directory with the verified clusters')
    parser.add_argument('finalized_table_path', help='path to an output file in which the final table will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        finalize_table(logger, args.putative_orthologs_path, args.verified_clusters_path, args.finalized_table_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
