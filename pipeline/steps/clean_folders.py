import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.general_utils import remove_path


def clean_folders(logger, job_input_path):
    with open(job_input_path) as fp:
        paths = [line.strip() for line in fp]

    for path in paths:
        remove_path(logger, path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path, help='path to file with the directories to clean')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, clean_folders, args.job_input_path)
