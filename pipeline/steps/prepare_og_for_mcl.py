
from sys import argv
import argparse
import logging
import itertools
import os
import sys
import traceback
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def add_hits_file_pair_to_mcl_inputs_files(logger, hits_file_path, gene_pair_to_og, ogs_texts_for_mcl):
    with open(hits_file_path, 'r') as fp:
        fp.readline()  # skip the header
        for line in fp:
            tokens = line.strip().split(',')
            if (tokens[0], tokens[1]) in gene_pair_to_og:
                og = gene_pair_to_og[(tokens[0], tokens[1])]
            elif (tokens[1], tokens[0]) in gene_pair_to_og:
                og = gene_pair_to_og[(tokens[1], tokens[0])]
            else:
                continue

            score = tokens[2]
            ogs_texts_for_mcl[og] += f'{tokens[0]}\t{tokens[1]}\t{score}\n'


def prepare_ogs_for_mcl(logger, hits_dir, putative_ogs_path, start, end, output_path):
    putative_ogs_df = pd.read_csv(putative_ogs_path, index_col=0)
    putative_ogs_df = putative_ogs_df.iloc[start:end]

    ogs_texts_for_mcl = {og_name: '' for og_name in putative_ogs_df.index}

    for hits_file_name in os.listdir(hits_dir):
        if not hits_file_name.endswith('.m8'):
            continue
        strain_1, strain_2 = os.path.splitext(hits_file_name)[0].split('_vs_')
        hits_file_path = os.path.join(hits_dir, hits_file_name)

        gene_pair_to_og = {}
        for og_name, og_row in putative_ogs_df.iterrows():
            if pd.isna(og_row[strain_1]) or pd.isna(og_row[strain_2]):
                continue

            strain_1_genes = og_row[strain_1].split(';')
            strain_2_genes = og_row[strain_2].split(';')

            if strain_1 != strain_2:
                for pair in itertools.product(strain_1_genes, strain_2_genes):
                    gene_pair_to_og[pair] = og_name
            else:
                for pair in itertools.combinations(strain_1_genes, 2):
                    gene_pair_to_og[pair] = og_name

        add_hits_file_pair_to_mcl_inputs_files(logger, hits_file_path, gene_pair_to_og, ogs_texts_for_mcl)

    for og_name, og_text in ogs_texts_for_mcl.items():
        mcl_file_path = os.path.join(output_path, f'{og_name}.mcl_input')
        if os.path.exists(mcl_file_path):
            return

        with open(mcl_file_path, 'w') as mcl_f:
            mcl_f.write(og_text)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('hits_dir', help='path to a dir with all hits')
    parser.add_argument('putative_ogs_path', help='path to a putative ogs path')
    parser.add_argument('start', help='first ortholog set to prepare', type=int)
    parser.add_argument('end', help='last ortholog set to prepare, exclusive', type=int)
    parser.add_argument('output_path', help='a folder in which the input files for mcl will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        prepare_ogs_for_mcl(logger, args.hits_dir, args.putative_ogs_path, args.start, args.end, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
