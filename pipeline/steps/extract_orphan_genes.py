#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 08:25:48 2023

@author: noa
"""
import json
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
from auxiliaries.plots_generator import generate_violinplot
from auxiliaries.logic_auxiliaries import get_all_genes_in_table


def get_all_genes_with_orthologs(orthologs_file):
    orthologs_df = pd.read_csv(orthologs_file)
    orthologs_df.drop(columns=['OG_name'], inplace=True)

    # Remove all rows that have only 1 strain in them (which are orphan genes that have multiple copies in the
    # same genome). Those rows will have not-nan values in 1 column
    all_clusters_no_orphans_df = orthologs_df[orthologs_df.count(axis=1) > 1]

    all_genes_with_orthologs = get_all_genes_in_table(all_clusters_no_orphans_df)
    return all_genes_with_orthologs


def extract_gene_names_from_fasta(orfs_file):
    with open(orfs_file) as ORFS:
        records = list(SeqIO.parse(ORFS, 'fasta'))
    gene_names = [record.name.strip() for record in records]
    return gene_names


def extract_orphan_proteins(logger, orfs_dir, orthologs_file, output_dir):
    if not os.path.exists(orfs_dir):
        logger.exception(f'ORFs dir does not exist in {orfs_dir}')
    if not os.path.exists(orthologs_file):
        logger.exception(f'Orthologs file does not exist in {orthologs_file}')
    if not os.path.exists(output_dir):
        logger.exception(f'Output path does not exist in {output_dir}')
    genes_with_orthologs = get_all_genes_with_orthologs(orthologs_file)
    number_of_orphans_per_file = {}
    for file_name in os.listdir(orfs_dir):
        org_name = os.path.splitext(file_name)[0]
        file_path = os.path.join(orfs_dir, file_name)
        if os.path.isfile(file_path):
            gene_names = set(extract_gene_names_from_fasta(file_path))
            orphans = list(gene_names.difference(genes_with_orthologs))
            number_of_orphans_per_file[org_name] = len(orphans)
            orphans_path = os.path.join(output_dir, f'{org_name}_orphans.txt')
            with open(orphans_path, 'w') as ORPH:
                ORPH.write('\n'.join(orphans))

    orphan_genes_count_per_genome_file_path = os.path.join(output_dir, 'orphan_genes_count.json')
    with open(orphan_genes_count_per_genome_file_path, 'w') as orphan_genes_count_per_genome_file:
        json.dump(number_of_orphans_per_file, orphan_genes_count_per_genome_file)

    orphan_genes_count_file_path = os.path.join(output_dir, 'orphan_genes_count.txt')
    with open(orphan_genes_count_file_path, 'w') as orphan_genes_count_file:
        orphan_genes_count_file.write('\n'.join([str(count) for count in number_of_orphans_per_file.values()]))

    generate_violinplot(orphan_genes_count_file_path, os.path.join(output_dir, 'orphan_genes_count.png'),
                        xlabel='Orphan genes count per genome', ylabel='Count')


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
        extract_orphan_proteins(logger, args.orfs_dir, args.orthologs_file, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')


if __name__ == '__main__':
    main()
