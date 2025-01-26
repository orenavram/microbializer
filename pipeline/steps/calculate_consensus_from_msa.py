import shutil
from sys import argv
import argparse
import logging
import subprocess
import os
import sys
import traceback
from Bio import AlignIO
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def calc_consensus(logger, og_path, output_dir):
    og_name = os.path.splitext(os.path.basename(og_path))[0]

    og_tmp_dir = os.path.join(output_dir, f'{og_name}_tmp')
    os.makedirs(og_tmp_dir, exist_ok=True)

    msa_stockholm_path = os.path.join(og_tmp_dir, f"{og_name}.sto")
    AlignIO.convert(og_path, "fasta", msa_stockholm_path, "stockholm")

    msa_db_path = os.path.join(og_tmp_dir, f"{og_name}_msaDb")
    profile_db_path = os.path.join(og_tmp_dir, f"{og_name}_profileDB")
    consensus_db_path = os.path.join(og_tmp_dir, f"{og_name}_consensusDb")
    consensus_faa_raw_path = os.path.join(og_tmp_dir, f"{og_name}_consensus.faa")

    cmds = [f"mmseqs convertmsa {msa_stockholm_path} {msa_db_path}",
            f"mmseqs msa2profile {msa_db_path} {profile_db_path} --match-mode 1",
            f"mmseqs profile2consensus {profile_db_path} {consensus_db_path}",
            f"mmseqs result2flat {consensus_db_path} {consensus_db_path} {consensus_db_path} {consensus_faa_raw_path}"]

    for cmd in cmds:
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True, check=True)

    record = SeqIO.parse(consensus_faa_raw_path, 'fasta').__next__()
    record.id = f"{og_name}_consensus"

    consensus_faa_path = os.path.join(output_dir, f"{og_name}.faa")
    SeqIO.write(record, consensus_faa_path, 'fasta')

    logger.info(f'Consensus calculation finished. Output written to {consensus_faa_path}')

    shutil.rmtree(og_tmp_dir)


def calc_consensus_from_all_ogs(logger, job_input_path, output_dir):
    with open(job_input_path, 'r') as f:
        ogs_paths = [line.strip() for line in f]

    for og_path in ogs_paths:
        calc_consensus(logger, og_path, output_dir)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', help='path to job input file with OG paths')
    parser.add_argument('output_dir', help='path to output dir')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        calc_consensus_from_all_ogs(logger, args.job_input_path, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
