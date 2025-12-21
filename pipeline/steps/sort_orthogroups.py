import argparse
from pathlib import Path
import sys
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step


def sort_orthogroups_df_and_rename_ogs(logger, orthogroups_file_path, orfs_coordinates_dir, sorted_orthogroups_file_path):
    def sort_genes_in_cell(cell, orfs_coordinates_dict):
        if pd.isna(cell) or cell == '':
            return cell

        genes = cell.split(';')
        genes_sorted = sorted(genes, key=lambda gene: orfs_coordinates_dict[gene])
        return ';'.join(genes_sorted)

    logger.info(f'Reading {orfs_coordinates_dir} into memory...')

    genome_name_to_orfs_coordinates = {}
    for orfs_coordinate_path in orfs_coordinates_dir.glob('*.csv'):
        genome_name = orfs_coordinate_path.stem
        orfs_coordinate_dict = pd.read_csv(orfs_coordinate_path, index_col=0)['coordinate'].to_dict()
        orfs_coordinate_dict[''] = float('inf')  # Add an empty string key with a value of inf to handle NaN values in the orthogroups DataFrame
        genome_name_to_orfs_coordinates[genome_name] = orfs_coordinate_dict

    logger.info(f'Finished reading {orfs_coordinates_dir} into memory.')

    # Sort the genes in each cell of the orthogroups DataFrame
    orthogroups_df = pd.read_csv(orthogroups_file_path, dtype=str)
    logger.info(f'Loaded orthogroups table from {orthogroups_file_path} into memory.')

    for col in orthogroups_df.columns[1:]:
        orthogroups_df[col] = orthogroups_df[col].apply(lambda cell: sort_genes_in_cell(cell, genome_name_to_orfs_coordinates[col]))
    logger.info('Sorted genes in each cell of the orthogroups table according to coordinates.')

    # Sort the table OGs by the first genome's coordinates, then by the second genome's coordinates, and so on.
    orthogroups_df = orthogroups_df.sort_values(
        by=list(orthogroups_df.columns[1:]),
        key=lambda col: col.map(lambda val: genome_name_to_orfs_coordinates[col.name][val.split(';')[0] if not pd.isna(val) else '']),
        ignore_index=True
    )
    logger.info('Sorted orthogroups table according to coordinates.')

    orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(orthogroups_df.index))]
    orthogroups_df.to_csv(sorted_orthogroups_file_path, index=False)
    logger.info(f'Wrote sorted orthogroups table according to coordinates to {sorted_orthogroups_file_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthogroups_file_path', type=Path)
    parser.add_argument('orfs_coordinates_dir', type=Path)
    parser.add_argument('sorted_orthogroups_file_path', type=Path)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, sort_orthogroups_df_and_rename_ogs, args.orthogroups_file_path, args.orfs_coordinates_dir,
             args.sorted_orthogroups_file_path)
