import argparse
from pathlib import Path
import sys
from time import sleep
import pandas as pd
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step


def load_all_ogs_hits(normalized_hits_dir, strain_to_gene_to_og):
    og_to_gene_pair_to_score = defaultdict(dict)

    for hits_file_path in normalized_hits_dir.glob('*.m8'):
        strain_1, strain_2 = hits_file_path.stem.split('_vs_')
        if strain_1 not in strain_to_gene_to_og or strain_2 not in strain_to_gene_to_og:  # it means that this hit file isn't relevant for the examined OGs
            continue

        strain_1_genes_to_og = strain_to_gene_to_og[strain_1]
        with open(hits_file_path) as f:
            f.readline()  # skip header
            for line in f:
                gene_1, gene_2, score = line.rstrip().split(',')
                # if gene_1 is in strain_1_genes_to_of and mapped to OG x, then gene_2 is also in OG x (becuase each putative OG is a conneceted component of the hits graph)
                if gene_1 in strain_1_genes_to_og:
                    og = strain_1_genes_to_og[gene_1]
                    og_to_gene_pair_to_score[og][(gene_1, gene_2)] = score

    return og_to_gene_pair_to_score


def prepare_ogs_for_mcl(logger, normalized_hits_dir, putative_ogs_path, job_input_path, output_path):
    with open(job_input_path, 'r') as f:
        ogs_numbers = [line.strip() for line in f]

    # Filter out OGs that already have an output file (this might happen if the job was restarted)
    ogs_numbers = [og_number for og_number in ogs_numbers if not (output_path / f'{og_number}.mcl_input').exists()]
    if not ogs_numbers:
        logger.info('All OGs already have an output file. Exiting...')
        return

    logger.info(f'Aggregating all genes from the specified {len(ogs_numbers)} putative OGs...')
    putative_ogs_df = pd.read_csv(putative_ogs_path, index_col=0)
    strain_to_gene_to_og = {}
    number_of_genes = 0
    for strain in putative_ogs_df.columns:
        strain_to_gene_to_og[strain] = {}
        for og_number in ogs_numbers:
            og_row = putative_ogs_df.loc[og_number]
            if pd.isna(og_row[strain]):
                continue
            for gene in og_row[strain].split(';'):
                strain_to_gene_to_og[strain][gene] = og_number
                number_of_genes += 1
    logger.info(f'All relevant genes were aggregated successfully. Number of relevant genes is {number_of_genes}.')

    logger.info('Loading relevant hits scores to memory...')
    og_to_gene_pair_to_score = load_all_ogs_hits(normalized_hits_dir, strain_to_gene_to_og)
    logger.info(f'All relevant hits were loaded successfully.')

    logger.info('Preparing input files for MCL...')
    for og_number in ogs_numbers:
        mcl_file_path = output_path / f'{og_number}.mcl_input'

        og_text_for_mcl = ''
        for (gene1, gene2), score in og_to_gene_pair_to_score[og_number].items():
            og_text_for_mcl += f'{gene1}\t{gene2}\t{score}\n'

        mcl_write_try_index = 1
        while not mcl_file_path.exists():
            try:
                with open(mcl_file_path, 'w') as mcl_f:
                    mcl_f.write(og_text_for_mcl)
                logger.info(f'Wrote {mcl_file_path}')
            except Exception as e:
                logger.error(f'Error writing {mcl_file_path} (try {mcl_write_try_index}): {e}')
                sleep(1)
                mcl_write_try_index += 1

    logger.info('Input files for MCL were written successfully.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('normalized_hits_dir', type=Path, help='path to a file with all hits')
    parser.add_argument('putative_ogs_path', type=Path, help='path to a putative ogs path')
    parser.add_argument('job_input_path', type=Path, help='')
    parser.add_argument('output_path', type=Path, help='a folder in which the input files for mcl will be written')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, prepare_ogs_for_mcl, args.normalized_hits_dir, args.putative_ogs_path, args.job_input_path,
             args.output_path)
