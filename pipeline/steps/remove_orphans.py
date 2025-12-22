import argparse
from pathlib import Path
import sys
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step


def remove_orphans(logger, orthogroups_file_path, output_orthogroups_file_path):
    orthogroups_df = pd.read_csv(orthogroups_file_path, index_col='OG_name', dtype=str)
    orthogroups_df = orthogroups_df[~(
            (orthogroups_df.count(axis=1) == 1) &
            ~(orthogroups_df.apply(lambda row: row.dropna().iloc[0].__contains__(';'), axis=1))
    )]
    orthogroups_df.reset_index(drop=True, inplace=True)
    orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(orthogroups_df.index))]
    orthogroups_df.set_index('OG_name', inplace=True)
    orthogroups_df.to_csv(output_orthogroups_file_path)
    logger.info(
        f'add_orphan_genes_to_ogs is False. Removed single orphan genes from {orthogroups_file_path} '
        f'and saved to {output_orthogroups_file_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthogroups_file_path', type=Path)
    parser.add_argument('output_orthogroups_file_path', type=Path)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, remove_orphans, args.orthogroups_file_path, args.output_orthogroups_file_path)
