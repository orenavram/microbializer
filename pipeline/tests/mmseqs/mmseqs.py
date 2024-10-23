import os
import subprocess
import sys
import time
from sys import argv
import argparse
import logging
import shutil
import pandas as pd
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries import consts
from auxiliaries.logic_auxiliaries import add_score_column_to_mmseqs_output

PROTEIN_FASTA_1 = 'Ecoli_042.02d_translated_orfs'
PROTEIN_FASTA_2 = 'Ecoli_536.02d_translated_orfs'


def main():
    m8_outfile = 'Ecoli_042_vs_Ecoli_536.m8'
    tmp_dir = Path('./tmp')

    cmd = f'mmseqs easy-rbh {PROTEIN_FASTA_1} {PROTEIN_FASTA_2} {m8_outfile} {tmp_dir} --format-output {consts.MMSEQS_OUTPUT_FORMAT} --threads 1'
    subprocess.run(cmd, shell=True)

    # Add 'score' column to mmseqs output
    df = pd.read_csv(m8_outfile, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER)
    add_score_column_to_mmseqs_output(df)

    df.to_csv(m8_outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()
