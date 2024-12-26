#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 08:25:48 2023

@author: noa
"""

from sys import argv
import argparse
import logging
import os
import sys
import pandas as pd
from Bio import SeqIO
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries.logic_auxiliaries import plot_genomes_histogram, flatten


def extract_gene_names_from_fasta(orfs_file):
    with open(orfs_file) as ORFS:
        gene_names = [record.name.strip() for record in SeqIO.parse(ORFS, 'fasta')]
    return gene_names


def extract_orphan_proteins(logger, orfs_file_path, orthogroups_file, output_dir):
    strain_name = os.path.splitext(os.path.basename(orfs_file_path))[0]

    orthogroups_df = pd.read_csv(orthogroups_file)
    orthogroups_df.drop(columns=['OG_name'], inplace=True)

    # Orphan ortogroups - orthogroups that contain genes from only one strain
    orphan_orthogroups = orthogroups_df[orthogroups_df.count(axis=1) == 1]
    orphan_orthogroups_of_strain = list(orphan_orthogroups[strain_name].dropna())
    genes_in_orphan_orthogroups_of_strain = flatten([orthogroup.split(';') for orthogroup in orphan_orthogroups_of_strain])

    # Orphan genes (not in any orthogroup)
    orthogroups_strain_column = orthogroups_df[strain_name].dropna()
    genes_in_orthogroups = flatten([value.split(';') for value in orthogroups_strain_column])
    all_strain_genes = set(extract_gene_names_from_fasta(orfs_file_path))
    orphans = list(all_strain_genes.difference(genes_in_orthogroups))

    orphans_path = os.path.join(output_dir, f'{strain_name}_orphans.txt')
    with open(orphans_path, 'w') as orphans_path_fp:
        orphans_path_fp.write('\n'.join(orphan_orthogroups_of_strain + orphans))

    orphans_count_path = os.path.join(output_dir, f'{strain_name}_orphans_stats.csv')
    orphans_stats = {
        'Orphan orthogroups count': len(orphan_orthogroups_of_strain),
        'Orphan single genes count': len(orphans),
        'Total orphans count': len(genes_in_orphan_orthogroups_of_strain) + len(orphans)
    }
    orphans_count_df = pd.DataFrame(orphans_stats, index=[strain_name])
    orphans_count_df.to_csv(orphans_count_path)


def extract_orphan_proteins_from_all_files(logger, job_input_path, orthogroups_file, output_dir):
    with open(job_input_path, 'r') as f:
        for line in f:
            orfs_file_path = line.strip()
            extract_orphan_proteins(logger, orfs_file_path, orthogroups_file, output_dir)


def main():
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', help='path to a file that contains the orfs path to extract orphan gene')
    parser.add_argument('orthogroups_file', help='path to the the orthologs table')
    parser.add_argument('output_dir', help='path to which the orphan proteins will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_orphan_proteins_from_all_files(logger, args.job_input_path, args.orthogroups_file, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)


if __name__ == '__main__':
    main()
