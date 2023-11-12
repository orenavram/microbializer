import pandas as pd
from sys import argv
import argparse
import os
import sys
import logging
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts


def get_orfs_lengths(query_vs_reference, orfs_statistics_dir):
    query_genome, subject_genome = os.path.splitext(os.path.basename(query_vs_reference))[0].split('_vs_')
    orfs_statistics_step_name = os.path.basename(orfs_statistics_dir)
    query_orfs_statistics_path = os.path.join(orfs_statistics_dir, f'{query_genome}.{orfs_statistics_step_name}')
    subject_orfs_statistics_path = os.path.join(orfs_statistics_dir, f'{subject_genome}.{orfs_statistics_step_name}')
    with open(query_orfs_statistics_path, 'r') as fp:
        query_orfs_statistics = json.load(fp)
    with open(subject_orfs_statistics_path, 'r') as fp:
        subject_orfs_statistics = json.load(fp)

    return query_orfs_statistics['lengths'], subject_orfs_statistics['lengths']


def filter_rbh_results(logger, query_vs_reference, output_path, precent_identity_cutoff, coverage_cutoff,
                       e_value_cutoff, orfs_statistics_dir, delimiter, names_delimiter):
    '''
    input:  path to file of blast results
            desired cutoff values
    output: file with filtered results
    '''
    logger.info(f'Filtering rbh results of {query_vs_reference}')

    query_orfs_lengths, subject_orfs_length = get_orfs_lengths(query_vs_reference, orfs_statistics_dir)

    # Here is the last time that '\t' is used as a delimiter!! from here and on, only ','
    df = pd.read_csv(query_vs_reference, header=None, sep='\t', names=consts.BLAST_OUTPUT_HEADERS)
    df['query_length'] = df['query'].apply(lambda query: query_orfs_lengths[query])
    df['subject_length'] = df['subject'].apply(lambda subject: subject_orfs_length[subject])
    df['query_coverage'] = df['query_end'] - (df['query_start']) / (df['query_length'] / 3)
    df['subject_coverage'] = df['alignment_length'] / (df['subject_length'] / 3)

    result = df[(df['identity_percent'] >= precent_identity_cutoff) & (df['evalue'] <= e_value_cutoff) &
                (df['query_coverage'] >= coverage_cutoff) & (df['subject_coverage'] >= coverage_cutoff)]
    columns_to_write = consts.BLAST_OUTPUT_HEADERS[:2] + consts.BLAST_OUTPUT_HEADERS[-1:]

    # e.g., ..../outputs/04_blast_filtered/Sflexneri_5_8401_vs_Ssonnei_Ss046.05_reciprocal_hits
    file_name = os.path.split(query_vs_reference)[-1]
    strain1_name, strain2_name = os.path.splitext(file_name)[0].split(names_delimiter)
    result.to_csv(output_path, sep=delimiter, index=False, header=[strain1_name, strain2_name, 'bitscore'],
                  columns=columns_to_write)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_result', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('--orfs_statistics_dir', help='path to directory of orfs statistics')
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=consts.CSV_DELIMITER)
    parser.add_argument('--names_delimiter', help='delimiter between the to species names', default='_vs_')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        filter_rbh_results(logger, args.blast_result, args.output_path, args.identity_cutoff, args.coverage_cutoff,
                           args.e_value_cutoff, args.orfs_statistics_dir, args.delimiter, args.names_delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
