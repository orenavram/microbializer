import argparse
from pathlib import Path
import sys
import pandas as pd
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.logic_utils import sort_orthogroups_by_columns


def merge_sub_orthogroups(logger, pseudo_orthogroups_file_path, sub_orthogroups_dir_path, output_path):
    temp_output_path = output_path.with_suffix('.tmp.csv')
    if temp_output_path.exists():
        df = pd.read_csv(temp_output_path, dtype=str, engine='pyarrow', dtype_backend='pyarrow')
        logger.info(f'Temporary merged sub-orthogroups file already exists at {temp_output_path}, so loaded it.')
    else:
        pseudo_orthogroups_df = pd.read_csv(pseudo_orthogroups_file_path, dtype=str, engine='pyarrow', dtype_backend='pyarrow')
        logger.info(f'Loaded pseudo-orthogroups file from {pseudo_orthogroups_file_path}')

        for sub_orthogroups_file_path in sorted(sub_orthogroups_dir_path.iterdir()):
            batch_id = sub_orthogroups_file_path.stem.split('_')[-1]
            logger.info(f'Merging orthogroups file from batch {batch_id} ({sub_orthogroups_file_path})...')

            sub_orthogroups_df = pd.read_csv(sub_orthogroups_file_path, dtype=str, engine='pyarrow', dtype_backend='pyarrow')
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
            logger.info(f'Merged orthogroups from batch {batch_id} into pseudo-orthogroups dataframe.')

        pseudo_orthogroups_df.to_csv(temp_output_path, index=False)
        logger.info(f'Merged sub-orthogroups file (unsorted) saved to {temp_output_path}')
        df = pseudo_orthogroups_df

    # Reorder columns
    sorted_cols = sorted(df.columns[1:])
    df = df[['OG_name'] + sorted_cols]

    # df = sort_orthogroups_by_columns(df)
    df = df.sort_values(by=list(df.columns[1:]), ignore_index=True)
    df.to_csv(output_path, index=False)
    logger.info(f'Merged sub-orthogroups file (sorted) saved to {output_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pseudo_orthogroups_file_path', type=Path)
    parser.add_argument('sub_orthogroups_dir_path', type=Path)
    parser.add_argument('output_path', type=Path)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, merge_sub_orthogroups, args.pseudo_orthogroups_file_path, args.sub_orthogroups_dir_path,
             args.output_path)
