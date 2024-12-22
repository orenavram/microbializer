import os
import subprocess
import sys
from sys import argv
import argparse
import logging
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger


def create_db(logger, protein_fasta_path, output_dir):
    strain_name = os.path.splitext(os.path.basename(protein_fasta_path))[0]
    db_path = os.path.join(output_dir, f'{strain_name}.db')
    if os.path.exists(db_path):
        return

    # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
    create_db_command = f'mmseqs createdb {protein_fasta_path} {db_path} --dbtype 1 -v 1'
    logger.info(f'Calling: {create_db_command}')
    subprocess.run(create_db_command, shell=True, check=True)

    logger.info(f"{db_path} was created successfully.")

    tmp_dir = os.path.join(output_dir, f'tmp_index_{strain_name}')
    create_index_command = f'mmseqs createindex {db_path} {tmp_dir} --threads 1 --search-type 1 --comp-bias-corr 0'
    logger.info(f'Calling: {create_index_command}')
    subprocess.run(create_index_command, shell=True, check=True)

    logger.info(f"Index for {db_path} was created successfully.")


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('protein_fasta_path', help='path to a protein fasta')
    parser.add_argument('output_dir', help='path to which the results will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        create_db(logger, args.protein_fasta_path, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
