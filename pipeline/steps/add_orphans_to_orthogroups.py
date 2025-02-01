from sys import argv
import argparse
import logging
import os
import sys
import pandas as pd
import re
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args


ORPHANS_FILENAME_GENOME_NAME_PATTERN = re.compile('(.+)_orphans.txt')


def finalize_table(logger, orthologs_table_path, finalized_table_path, orphan_genes_dir):
    orthogroups_df = pd.read_csv(orthologs_table_path)
    logger.info(f'Read {orthologs_table_path} into memory. There are {len(orthogroups_df.index)} orthogroups.')

    logger.info(f'Starting to aggregate orphan genes from {orphan_genes_dir}')
    orphan_genes = []
    for filename in os.listdir(orphan_genes_dir):
        strain_match_object = ORPHANS_FILENAME_GENOME_NAME_PATTERN.match(filename)
        if not strain_match_object:
            continue
        strain = strain_match_object.group(1)
        logger.info(f'Aggregating orphan genes from {strain}...')
        with open(os.path.join(orphan_genes_dir, filename)) as orphan_genes_file:
            # add only single orphans and not orthogroup orphans since those already are in the orthogroups table.
            orphan_genes_to_add = [line.strip() for line in orphan_genes_file if line and ';' not in line]

        orphan_genes.extend([(strain, gene) for gene in orphan_genes_to_add])
        logger.info(f'Finished aggregating orphan genes from {strain}')

    logger.info(f"Found {len(orphan_genes)} orphan genes to add as new OGs")
    orphan_clusters_df = pd.DataFrame([pd.Series({'OG_name': '', strain: gene}) for strain, gene in orphan_genes])
    orthogroups_df = pd.concat([orthogroups_df, orphan_clusters_df], ignore_index=True)
    logger.info(f'Finished adding orphan genes as OGs. OG table now contains {len(orthogroups_df.index)} groups.')

    # Sort the rows by all columns except the first one (OG_name) to keep consistent output
    orthogroups_df = orthogroups_df.sort_values(by=list(orthogroups_df.columns[1:])).reset_index(drop=True)
    orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(orthogroups_df.index))]
    orthogroups_df.to_csv(finalized_table_path, index=False)
    logger.info(f'Wrote finalized sorted orthogroups table to {finalized_table_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', help='path to the orthologs table (input)')
    parser.add_argument('finalized_table_path', help='path to the output orthologs table')
    parser.add_argument('--orphan_genes_dir', help='path to a directory with the orphan genes')
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        finalize_table(logger, args.orthologs_table_path, args.finalized_table_path, args.orphan_genes_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
