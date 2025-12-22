import argparse
from pathlib import Path
import sys
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.logic_utils import split_ogs_to_jobs_inputs_files_by_og_sizes


def split_ogs_by_size(logger, orthogroups_file_path, jobs_inputs_dir, max_parallel_jobs):
    orthogroups_df = pd.read_csv(orthogroups_file_path, dtype=str)
    logger.info(f'Read orthogroups file from {orthogroups_file_path} into memory')

    split_ogs_to_jobs_inputs_files_by_og_sizes(orthogroups_df, jobs_inputs_dir, max_parallel_jobs)
    logger.info(f'Split orthogroups into {max_parallel_jobs} jobs input files in {jobs_inputs_dir} based on OG sizes')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthogroups_file_path', type=Path)
    parser.add_argument('jobs_inputs_dir', type=Path)
    parser.add_argument('max_parallel_jobs', type=int)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, split_ogs_by_size, args.orthogroups_file_path, args.jobs_inputs_dir, args.max_parallel_jobs)
