from sys import argv
import argparse
from pathlib import Path
import sys
import traceback
import pandas as pd
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args


def merge_sub_orthogroups(logger, pseudo_orthogroups_file_path, sub_orthogroups_dir_path, output_path):
    pseudo_orthogroups_df = pd.read_csv(pseudo_orthogroups_file_path)

    for sub_orthogroups_file_path in sub_orthogroups_dir_path.iterdir():
        batch_id = sub_orthogroups_file_path.stem.split('_')[-1]

        sub_orthogroups_df = pd.read_csv(sub_orthogroups_file_path)
        sub_orthogroups_df.drop(columns=['OG_name'], inplace=True)
        sub_orthogroups_df.set_index('representative_gene', inplace=True)

        # Fill rows of pseudo_orthogroups_df where pseudo_genome_{batch_id} contains only 1 gene
        pseudo_orthogroups_df = pseudo_orthogroups_df.merge(
            sub_orthogroups_df, how='left', left_on=f'pseudo_genome_{batch_id}', right_on='representative_gene')

        # Fill rows of pseudo_orthogroups_df where pseudo_genome_{batch_id} contains multiple genes
        pseudo_orthogroups_with_paralogs_in_pseudo_genome = pseudo_orthogroups_df[
            pseudo_orthogroups_df[f'pseudo_genome_{batch_id}'].str.contains(';', na=False)]
        for i, row in pseudo_orthogroups_with_paralogs_in_pseudo_genome.iterrows():
            pseudo_genes = row[f'pseudo_genome_{batch_id}'].split(';')
            row_sub_orthogroups = pd.concat([sub_orthogroups_df.loc[pseudo_gene] for pseudo_gene in pseudo_genes],
                                            axis=1)
            strain_to_genes = dict(
                row_sub_orthogroups.apply(lambda row: ';'.join(sorted(row.dropna().astype(str))), axis=1))
            for strain, genes in strain_to_genes.items():
                pseudo_orthogroups_df.at[i, strain] = genes if genes else np.nan

        pseudo_orthogroups_df.drop(columns=[f'pseudo_genome_{batch_id}'], inplace=True)

    pseudo_orthogroups_df = pseudo_orthogroups_df.sort_values(
        by=list(pseudo_orthogroups_df.columns[1:])).reset_index(drop=True)
    pseudo_orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(pseudo_orthogroups_df.index))]
    pseudo_orthogroups_df.to_csv(output_path, index=False)
    logger.info(f'Merged sub-orthogroups file saved to {output_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('pseudo_orthogroups_file_path', type=Path)
    parser.add_argument('sub_orthogroups_dir_path', type=Path)
    parser.add_argument('output_path', type=Path)
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        merge_sub_orthogroups(logger, args.pseudo_orthogroups_file_path, args.sub_orthogroups_dir_path, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
