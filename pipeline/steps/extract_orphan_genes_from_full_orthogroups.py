#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 08:25:48 2023

@author: noa
"""

from sys import argv
import argparse
from pathlib import Path
import os
import sys
import pandas as pd
from Bio import SeqIO
import traceback

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger, add_default_step_args
from auxiliaries.logic_auxiliaries import plot_genomes_histogram, flatten


def extract_gene_names_from_fasta(proteins_file_path):
    with open(proteins_file_path) as proteins_file_fp:
        gene_names = [record.name.strip() for record in SeqIO.parse(proteins_file_fp, 'fasta')]
    return gene_names


def extract_orphan_proteins(logger, strain_name, orphan_orthogroups, output_dir):
    orphan_orthogroups_of_strain = orphan_orthogroups[strain_name].dropna()

    orphan_orthogroups_with_paralogs = orphan_orthogroups_of_strain[orphan_orthogroups_of_strain.str.contains(';')]
    genes_in_orphan_orthogroups_with_paralogs = flatten(
        [orthogroup.split(';') for orthogroup in orphan_orthogroups_with_paralogs])

    orphan_genes = orphan_orthogroups_of_strain[~orphan_orthogroups_of_strain.str.contains(';')]

    orphans_path = os.path.join(output_dir, f'{strain_name}_orphans.txt')
    with open(orphans_path, 'w') as orphans_path_fp:
        orphans_path_fp.write('\n'.join(list(orphan_orthogroups_with_paralogs) + list(orphan_genes)))
    logger.info(f'Orphan genes of {strain_name} were written to {orphans_path}')

    orphans_count_path = os.path.join(output_dir, f'{strain_name}_orphans_stats.csv')
    orphans_stats = {
        'Orphan orthogroups count': len(orphan_orthogroups_with_paralogs),
        'Orphan single genes count': len(orphan_genes),
        'Total orphans count': len(genes_in_orphan_orthogroups_with_paralogs) + len(orphan_genes)
    }
    orphans_count_df = pd.DataFrame(orphans_stats, index=[strain_name])
    orphans_count_df.to_csv(orphans_count_path)
    logger.info(f'Orphan genes statistics of {strain_name} were written to {orphans_count_path}')


def extract_orphan_proteins_of_all_strains(logger, job_input_path, orthogroups_file, output_dir):
    # Here we start from orthogroups_df that already contains orthogroups for all orphan genes.
    orthogroups_df = pd.read_csv(orthogroups_file)
    orthogroups_df.drop(columns=['OG_name'], inplace=True)
    orphan_orthogroups = orthogroups_df[orthogroups_df.count(axis=1) == 1]

    with open(job_input_path, 'r') as f:
        for line in f:
            strain_name = line.strip()
            extract_orphan_proteins(logger, strain_name, orphan_orthogroups, output_dir)


def main():
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path, help='path to a file that contains the strain names to extract orphan genes')
    parser.add_argument('orthogroups_file', type=Path, help='path to the the orthologs table')
    parser.add_argument('output_dir', type=Path, help='path to which the orphan proteins will be written')
    add_default_step_args(parser)
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose)

    logger.info(script_run_message)
    try:
        extract_orphan_proteins_of_all_strains(logger, args.job_input_path, args.orthogroups_file, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)


if __name__ == '__main__':
    main()
