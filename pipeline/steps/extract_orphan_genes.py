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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries.logic_auxiliaries import get_all_genes_in_table, plot_genomes_histogram, flatten


def extract_gene_names_from_fasta(orfs_file):
    with open(orfs_file) as ORFS:
        records = list(SeqIO.parse(ORFS, 'fasta'))
    gene_names = [record.name.strip() for record in records]
    return gene_names


def extract_orphan_proteins(logger, orfs_dir, orthologs_file, output_dir):
    orthologs_df = pd.read_csv(orthologs_file)
    orthologs_df.drop(columns=['OG_name'], inplace=True)

    orthogroups_with_one_strain = orthologs_df[orthologs_df.count(axis=1) == 1]

    all_genes_in_orthogroups = get_all_genes_in_table(orthologs_df)

    number_of_orphans_per_file = {}
    for orfs_file_name in os.listdir(orfs_dir):
        species_name = os.path.splitext(orfs_file_name)[0]
        orfs_file_path = os.path.join(orfs_dir, orfs_file_name)
        if os.path.isfile(orfs_file_path):
            gene_names = set(extract_gene_names_from_fasta(orfs_file_path))
            orphan_orthogroups = list(orthogroups_with_one_strain[species_name].dropna())
            orphans = list(gene_names.difference(all_genes_in_orthogroups))
            number_of_orphans_per_file[species_name] = len(flatten([orthogroup.split(';') for orthogroup in orphan_orthogroups])) + len(orphans)
            orphans_path = os.path.join(output_dir, f'{species_name}_orphans.txt')
            with open(orphans_path, 'w') as orphans_path_fp:
                orphans_path_fp.write('\n'.join(orphan_orthogroups))
                orphans_path_fp.write('\n'.join(orphans))

    plot_genomes_histogram(number_of_orphans_per_file, output_dir, 'orphan_genes_count', 'Orphan genes count',
                           'Orphan genes count per Genome')


def main():
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('orfs_dir', help='path to the ORFs dir')
    parser.add_argument('orthologs_file', help='path to the the orthologs table')
    parser.add_argument('output_dir', help='path to which the orphan proteins will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        if not os.path.exists(args.orfs_dir):
            logger.exception(f'ORFs dir does not exist in {args.orfs_dir}')
        if not os.path.exists(args.orthologs_file):
            logger.exception(f'Orthologs file does not exist in {args.orthologs_file}')
        if not os.path.exists(args.output_dir):
            logger.exception(f'Output path does not exist in {args.output_dir}')

        extract_orphan_proteins(logger, args.orfs_dir, args.orthologs_file, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')


if __name__ == '__main__':
    main()
