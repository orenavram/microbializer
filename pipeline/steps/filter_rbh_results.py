import numpy as np
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


def filter_rbh_results(logger, query_vs_reference, output_path, scores_statistics_dir, precent_identity_cutoff,
                       coverage_cutoff, e_value_cutoff, delimiter, names_delimiter):
    '''
    input:  path to file of blast results
            desired cutoff values
    output: file with filtered results
    '''
    logger.info(f'Filtering rbh results of {query_vs_reference}')

    # Here is the last time that '\t' is used as a delimiter!! from here and on, only ','
    df = pd.read_csv(query_vs_reference, sep='\t')

    result = df[(df['fident'] >= precent_identity_cutoff) & (df['evalue'] <= e_value_cutoff) &
                (df['qcov'] >= coverage_cutoff) & (df['tcov'] >= coverage_cutoff)]

    # e.g., ..../outputs/04_blast_filtered/Sflexneri_5_8401_vs_Ssonnei_Ss046.05_reciprocal_hits
    query_vs_reference_file_name = os.path.splitext(os.path.basename(query_vs_reference))[0]
    strain1_name, strain2_name = query_vs_reference_file_name.split(names_delimiter)
    result.to_csv(output_path, sep=delimiter, index=False, header=[strain1_name, strain2_name, 'score'],
                  columns=['query', 'target', 'score'])

    score_stats_file = os.path.join(scores_statistics_dir, f'{query_vs_reference_file_name}.stats')
    statistics = {'mean': np.mean(result['score']), 'sum': np.sum(result['score']),
                  'number of records': len(result['score'])}
    with open(score_stats_file, 'w') as fp:
        json.dump(statistics, fp)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_result', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('scores_statistics_dir', help='path to output dir of score statistics')
    parser.add_argument('--identity_cutoff', type=float)
    parser.add_argument('--coverage_cutoff', type=float)
    parser.add_argument('--e_value_cutoff', type=float)
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=consts.CSV_DELIMITER)
    parser.add_argument('--names_delimiter', help='delimiter between the to species names', default='_vs_')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        filter_rbh_results(logger, args.blast_result, args.output_path, args.scores_statistics_dir, args.identity_cutoff,
                           args.coverage_cutoff, args.e_value_cutoff, args.delimiter, args.names_delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
