from sys import argv
import argparse
from pathlib import Path
import sys
import pandas as pd
from collections import defaultdict
import traceback

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args
from auxiliaries.logic_auxiliaries import get_strain_name


def get_verified_clusters_set(verified_clusters_path):
    return set([file.stem for file in verified_clusters_path.glob('*verified_cluster')])


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
    for cluster_file_path in verified_clusters_path.glob('*split_cluster'):
        with open(cluster_file_path, 'r') as new_cluster:
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
    parser.add_argument('putative_orthologs_path', type=Path, help='path to a file with the putative orthologs sets')
    parser.add_argument('verified_clusters_path', type=Path, help='path to a directory with the verified clusters')
    parser.add_argument('finalized_table_path', type=Path, help='path to an output file in which the final table will be written')
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        finalize_table(logger, args.putative_orthologs_path, args.verified_clusters_path, args.finalized_table_path)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
