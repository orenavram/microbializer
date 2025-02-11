import argparse
import sys
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import add_default_step_args, run_step


def finalize_table(logger, orthologs_table_path, finalized_table_path, orphan_genes_dir):
    orthogroups_df = pd.read_csv(orthologs_table_path)
    logger.info(f'Read {orthologs_table_path} into memory. There are {len(orthogroups_df.index)} orthogroups.')

    logger.info(f'Starting to aggregate orphan genes from {orphan_genes_dir}')
    orphan_genes = []
    for orphan_genes_file_path in orphan_genes_dir.glob('*_orphans.txt'):
        strain = str(orphan_genes_file_path.stem).replace('_orphans', '')
        logger.info(f'Aggregating orphan genes from {strain}...')
        with open(orphan_genes_file_path) as orphan_genes_file:
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
    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', type=Path, help='path to the orthologs table (input)')
    parser.add_argument('finalized_table_path', type=Path, help='path to the output orthologs table')
    parser.add_argument('--orphan_genes_dir', type=Path, help='path to a directory with the orphan genes')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, finalize_table, args.orthologs_table_path, args.finalized_table_path, args.orphan_genes_dir)
