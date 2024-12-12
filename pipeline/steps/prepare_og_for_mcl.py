
from sys import argv
import argparse
import logging
from itertools import combinations
import os
import sys
import traceback
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts


def prepare_og_for_mcl(logger, hits_dir, putative_og_path, output_path):
    putative_og_name = os.path.basename(putative_og_path).split('.')[0]

    mcl_file_path = os.path.join(output_path, putative_og_name + '.mcl_input')
    if os.path.exists(mcl_file_path):
        return

    putative_og_df = pd.read_csv(putative_og_path)
    putative_og_df.drop(columns=['OG_name'], inplace=True)

    gene_to_strain = {}
    for col in putative_og_df.columns:
        col_genes = putative_og_df[col].values[0]
        if pd.notna(col_genes):
            for gene in col_genes.split(';'):
                gene_to_strain[gene] = col

    text_to_mcl_file = ''
    for gene_1, gene_2 in set(combinations(gene_to_strain.keys(), 2)):
        strain_1 = gene_to_strain[gene_1]
        strain_2 = gene_to_strain[gene_2]

        hits_path = os.path.join(hits_dir, f'{strain_1}_vs_{strain_2}.m8')
        if not os.path.exists(hits_path):
            hits_path = os.path.join(hits_dir, f'{strain_2}_vs_{strain_1}.m8')
            if not os.path.exists(hits_path):
                logger.warning(f'No hits file for {strain_1} vs {strain_2}')
                continue

        hits_df = pd.read_csv(hits_path)
        relevant_hit = hits_df[(hits_df[strain_1] == gene_1) & (hits_df[strain_2] == gene_2)]
        if relevant_hit.empty:
            # genes pair wasn't a hit
            continue

        score = relevant_hit['score'].values[0]
        text_to_mcl_file += f'{gene_1}\t{gene_2}\t{score}\n'

    with open(mcl_file_path, 'w') as mcl_f:
        mcl_f.write(text_to_mcl_file)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('hits_dir', help='path to a dir with all hits')
    parser.add_argument('putative_og_path', help='path to a putative og file')
    parser.add_argument('output_path', help='a folder in which the input files for mcl will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        prepare_og_for_mcl(logger, args.hits_dir, args.putative_og_path, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
