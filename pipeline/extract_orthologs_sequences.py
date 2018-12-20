# -*- coding: utf-8 -*-
"""
Created on Tue May 22 14:52:18 2018

@author: shirportugez

This script is a module that extract sequences of given ortholog groups.
Input:
    (1) a path to orthologs table.
    (2) directory for sequence files (each bacteria in the og table should have one sequence file with the name of the bacteria that should contain all genes of that bacteria).
    (3) a directory for output sequence files - if the folder doesn't exists, the script will create it.
    (4) optional [fasta by default]: input and output sequence file format (fasta or genebank)
Output:
    directory containing fasta files for each ortholog group.
"""

import os

def get_orthologs_group_sequences(orthologs_group, sequences_dir):
    pass


def extract_orthologs_sequences(orthologs_table_path, sequences_dir, output_path, number_of_first_group, delimiter):
    i = number_of_first_group
    with open(orthologs_table_path) as f:
        for line in f:
            orthologs_group = line.rstrip().split(delimiter)
            orthologs_sequences = get_orthologs_group_sequences(orthologs_group, sequences_dir)
            with open(os.path.join(output_path, i+'_orthologs_group.fasta'), 'w'):
                f.write(orthologs_sequences)
            i += 1


if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    logger.info(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', help='path to an orthologs table (aka verified clusters)')
    parser.add_argument('sequences_dir', help='path to a directory with the bacterial gene sequences (aka ORFs)')
    parser.add_argument('output_path', help='path to an output directory (aka orthologs sets sequences)')
    parser.add_argument('--number_of_first_group', help='number that represents the "arbirary" name that the first '
                        'orthologs group will get (the second will get number_of_first_group + 1, etc...)', default = '0')
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=',')
    args = parser.parse_args()

    extract_orthologs_sequences(args.orthologs_table_path, args.sequences_dir, args.output_path, args.number_of_first_group, args.delimiter)
