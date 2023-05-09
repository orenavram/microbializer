import pandas as pd
from sys import argv
import argparse
import os
import sys
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def filter_rbh_results(logger, query_vs_reference, output_path, precent_identity_cutoff,
                       e_value_cutoff, delimiter, names_delimiter):
    '''
    input:  path to file of blast results
            desired cutoff values
    output: file with filtered results
    '''
    logger.info(f'Filtering rbh results of {query_vs_reference}')

    header = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    # Here is the last time that '\t' is used as a delimiter!! from here and on, only ','
    df = pd.read_csv(query_vs_reference, header=None, sep='\t', names=header.split())
    result = df[(df.pident >= precent_identity_cutoff) & (df.evalue <= e_value_cutoff)]
    splitted_header = header.split()
    columns_to_write = splitted_header[:2] + splitted_header[-1:]

    # e.g., ..../outputs/04_blast_filtered/Sflexneri_5_8401_vs_Ssonnei_Ss046.05_reciprocal_hits
    file_name = os.path.split(query_vs_reference)[-1]
    strain1_name, strain2_name = os.path.splitext(file_name)[0].split(names_delimiter)
    result.to_csv(output_path, sep=delimiter, index=False, header=[strain1_name, strain2_name, 'bitscore'],
                  columns=columns_to_write)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('blast_result', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('--identity_cutoff', help='path to translated sequences', type=float)
    parser.add_argument('--e_value_cutoff', help='path to translated sequences', type=float)
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=',')
    parser.add_argument('--names_delimiter', help='delimiter between the to species names', default='_vs_')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        filter_rbh_results(logger, args.blast_result, args.output_path, args.identity_cutoff,
                           args.e_value_cutoff, args.delimiter, args.names_delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
