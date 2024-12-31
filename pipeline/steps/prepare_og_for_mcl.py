
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
from auxiliaries.logic_auxiliaries import flatten


def load_all_ogs_hits(all_reciprocal_hits_path, all_relevant_genes):
    gene_pair_to_score = {}
    with open(all_reciprocal_hits_path) as f:
        for line in f:
            line_tokens = line.rstrip().split(',')
            if 'score' in line:
                # new reciprocal hits file starts
                continue
            gene_1, gene_2 = line_tokens[:2]
            if gene_1 in all_relevant_genes:  # gene_1 in all_relevant_genes <==> gene_2 in all_relevant_genes
                score = line_tokens[2]
                gene_pair_to_score[(gene_1, gene_2)] = score

    return gene_pair_to_score


def prepare_ogs_for_mcl(logger, all_reciprocal_hits_path, putative_ogs_path, job_input_path, output_path):
    putative_ogs_df = pd.read_csv(putative_ogs_path, index_col=0)
    with open(job_input_path, 'r') as f:
        ogs_numbers = [f'OG_{num}' for num in f.read().splitlines()]

    logger.info(f'Aggregating all genes from the specified {len(ogs_numbers)} putative OGs...')
    all_relevant_genes = set()
    og_to_genes = {}
    for og_number in ogs_numbers:
        og_row = putative_ogs_df.loc[og_number]
        og_genes = flatten([strain_genes.split(';') for strain_genes in og_row if not pd.isna(strain_genes)])
        og_to_genes[og_number] = og_genes
        all_relevant_genes.update(og_genes)
    logger.info(f'All relevant genes were aggregated successfully. Number of relevant genes is {len(all_relevant_genes)}.')

    logger.info('Loading relevant hits scores to memory...')
    gene_pair_to_score = load_all_ogs_hits(all_reciprocal_hits_path, all_relevant_genes)
    logger.info(f'All relevant hits were loaded successfully. Number of relevant hits (gene pairs) is {len(gene_pair_to_score)}.')

    logger.info('Preparing input files for MCL...')
    for og_number, genes in og_to_genes.items():
        mcl_file_path = os.path.join(output_path, og_number + '.mcl_input')
        if os.path.exists(mcl_file_path):
            continue

        og_text_for_mcl = ''
        for gene1, gene2 in itertools.combinations(genes, 2):
            score = gene_pair_to_score.get((gene1, gene2), gene_pair_to_score.get((gene2, gene1)))
            if score:
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
