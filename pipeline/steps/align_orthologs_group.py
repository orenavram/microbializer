from sys import argv
import argparse
import logging
import subprocess
import os
import sys
import traceback
import shutil
from Bio import SeqIO
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def reconstruct_msa(logger, sequences_file_path, mode, output_file_path):
    # cmd = f'mafft-linsi --maxiterate {maxiterate} --localpair {sequences_file_path} > {output_file_path}'

    # --auto Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i and FFT-NS-2, according to data size.
    # --amino/--nuc tells mafft that's an amino acid/nucleotide (respectively) msa. If you let it decide by itself, it
    # might wrong on small data sets as they might look like dna but they are NOT! e.g.,
    # [orenavr2@powerlogin-be2 test]$ cat /groups/pupko/orenavr2/igomeProfilingPipeline/experiments/test/analysis/motif_inference/17b_03/unaligned_sequences/17b_03_clusterRank_215_uniqueMembers_2_clusterSize_252.81.faa
    # >seq_235_lib_12_len_12_counts_126.40626975097965
    # CNTDVACAAPGN
    # >seq_1112_lib_C8C_len_10_counts_126.40626975097965
    # CTTACAPVNC
    start_time = time.time()
    records = list(SeqIO.parse(sequences_file_path, 'fasta'))
    parse_time = time.time() - start_time
    logger.info(f'Biopython parse time is {parse_time}')
    if len(records) == 1:
        logger.info(f'Only one sequence in {sequences_file_path}. Copying it to {output_file_path}')
        shutil.copy(sequences_file_path, output_file_path)
    else:
        start_time = time.time()
        cmd = f'mafft --auto --{mode} {sequences_file_path} > {output_file_path}'
        logger.info(f'Starting MAFFT. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)
        mafft_time = time.time() - start_time
        logger.info(f'mafft time is {mafft_time}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('sequences_file_path', help='path to a file with unaligned sequences')
    parser.add_argument('output_file_path', help='path to a file in which the aligned sequences will be written')
    parser.add_argument('--type', choices=['amino', 'nuc'], default='amino',
                        help="amino/nuc tells mafft that's an amino acid/nucleotide (respectively) msa")
    # parser.add_argument('--maxiterate', help='number for MAFFT maxiterate parameter', default=1000, type=int)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        reconstruct_msa(logger, args.sequences_file_path, args.type, args.output_file_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
