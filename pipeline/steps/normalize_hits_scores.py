import pandas as pd
import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.general_utils import str_to_bool


def normalize_hits_scores(logger, blast_result, output_path, scores_normalize_coefficient, use_parquet):
    query_vs_reference_file_name = blast_result.stem
    strain1_name, strain2_name = query_vs_reference_file_name.split('_vs_')

    if use_parquet:
        df = pd.read_parquet(blast_result)
    else:
        df = pd.read_csv(blast_result)

    df['score'] = (df['score'] / scores_normalize_coefficient).round(2)
    df.to_csv(output_path, index=False, header=[strain1_name, strain2_name, 'score'])
    logger.info(f'Normalized scores of {blast_result} written to {output_path}')


def normalize_hits_scores_of_all_files(logger, job_input_path, output_dir, use_parquet):
    with open(job_input_path, 'r') as f:
        for line in f:
            hits_path, scores_normalize_coefficient = line.strip().split()
            hits_path = Path(hits_path)
            output_path = output_dir / f"{hits_path.stem}.m8"
            normalize_hits_scores(logger, hits_path, output_path, float(scores_normalize_coefficient), use_parquet)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir', type=Path, help='path to output dir')
    parser.add_argument('--use_parquet', type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, normalize_hits_scores_of_all_files, args.job_input_path, args.output_dir, args.use_parquet)
