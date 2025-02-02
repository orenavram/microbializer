import os
import subprocess
import sys
from sys import argv
import argparse
from pathlib import Path
import traceback

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, add_default_step_args


def create_db(logger, genome_name, proteomes_dir, output_dir):
    protein_fasta_path = os.path.join(proteomes_dir, f'{genome_name}.faa')
    db_path = os.path.join(output_dir, f'{genome_name}.db')
    if os.path.exists(db_path):
        return

    # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
    create_db_command = f'mmseqs createdb {protein_fasta_path} {db_path} --dbtype 1 -v 1'
    logger.info(f'Calling: {create_db_command}')
    subprocess.run(create_db_command, shell=True, check=True)

    logger.info(f"{db_path} was created successfully.")

    # tmp_dir = os.path.join(output_dir, f'tmp_index_{genome_name}')
    # create_index_command = f'mmseqs createindex {db_path} {tmp_dir} --threads 1 --search-type 1 --comp-bias-corr 0'
    # logger.info(f'Calling: {create_index_command}')
    # subprocess.run(create_index_command, shell=True, check=True)

    # logger.info(f"Index for {db_path} was created successfully.")


def create_dbs(logger, job_input_path, proteomes_dir, output_dir):
    with open(job_input_path) as fp:
        for line in fp:
            genome_name = line.strip()
            create_db(logger, genome_name, proteomes_dir, output_dir)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path, help='path to a file that contains the genome names to create dbs for')
    parser.add_argument('proteomes_dir', type=Path, help='path to dir of proteomes')
    parser.add_argument('output_dir', type=Path, help='path to which the results will be written')
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        create_dbs(logger, args.job_input_path, args.proteomes_dir, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
