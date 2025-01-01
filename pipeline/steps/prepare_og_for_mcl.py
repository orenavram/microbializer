
from sys import argv
import argparse
import logging
import os
import sys
import traceback
import pandas as pd
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries.logic_auxiliaries import flatten


def load_all_ogs_hits(all_reciprocal_hits_path, gene_to_og):
    og_to_gene_pair_to_score = defaultdict(dict)
    with open(all_reciprocal_hits_path) as f:
        for line in f:
            line_tokens = line.rstrip().split(',')
            if 'score' in line:
                # new reciprocal hits file starts
                continue
            gene_1, gene_2 = line_tokens[:2]
            if gene_1 in gene_to_og:  # if gene_1 is in gene_to_og and mapped to OG x, then gene_2 is also in OG x
                                      # (becuase each putative OG is a conneceted component of the graph described in all_reciprocal_hits_path)
                og = gene_to_og[gene_1]
                score = line_tokens[2]
                og_to_gene_pair_to_score[og][(gene_1, gene_2)] = score

    return og_to_gene_pair_to_score


def prepare_ogs_for_mcl(logger, all_reciprocal_hits_path, putative_ogs_path, job_input_path, output_path):
    putative_ogs_df = pd.read_csv(putative_ogs_path, index_col=0)
    with open(job_input_path, 'r') as f:
        ogs_numbers = [f'OG_{num}' for num in f.read().splitlines()]

    logger.info(f'Aggregating all genes from the specified {len(ogs_numbers)} putative OGs...')
    # all_relevant_genes = set()
    # og_to_genes = {}
    gene_to_og = {}
    for og_number in ogs_numbers:
        og_row = putative_ogs_df.loc[og_number]
        og_genes = flatten([strain_genes.split(';') for strain_genes in og_row if not pd.isna(strain_genes)])
        # og_to_genes[og_number] = og_genes
        # all_relevant_genes.update(og_genes)
        for gene in og_genes:
            gene_to_og[gene] = og_number
    logger.info(f'All relevant genes were aggregated successfully. Number of relevant genes is {len(gene_to_og)}.')

    logger.info('Loading relevant hits scores to memory...')
    og_to_gene_pair_to_score = load_all_ogs_hits(all_reciprocal_hits_path, gene_to_og)
    logger.info(f'All relevant hits were loaded successfully.')

    logger.info('Preparing input files for MCL...')
    for og_number in ogs_numbers:
        mcl_file_path = os.path.join(output_path, og_number + '.mcl_input')
        if os.path.exists(mcl_file_path):
            continue

        og_text_for_mcl = ''
        for (gene1, gene2), score in og_to_gene_pair_to_score[og_number].items():
            og_text_for_mcl += f'{gene1}\t{gene2}\t{score}\n'

        with open(mcl_file_path, 'w') as mcl_f:
            mcl_f.write(og_text_for_mcl)

    logger.info('Input files for MCL were written successfully.')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('all_reciprocal_hits_path', help='path to a file with all hits')
    parser.add_argument('putative_ogs_path', help='path to a putative ogs path')
    parser.add_argument('job_input_path', help='')
    parser.add_argument('output_path', help='a folder in which the input files for mcl will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        prepare_ogs_for_mcl(logger, args.all_reciprocal_hits_path, args.putative_ogs_path, args.job_input_path, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
