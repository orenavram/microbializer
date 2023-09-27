import sys
import os
from sys import argv
import argparse
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import load_header2sequences_dict, get_job_logger


def induce_sequence(logger, aa_seq, dna_seq):
    result = ''
    dna_i = 0
    for aa_i in range(len(aa_seq)):
        if aa_seq[aa_i] == '-':
            result += '-' * 3
        else:
            result += dna_seq[dna_i:dna_i + 3]
            dna_i += 3

    # TODO: remove this checkup
    if len(aa_seq) * 3 != len(result):
        logger.error('$' * 80)
        logger.error('len(aa_seq)*3 != len(result)')
        logger.error(f'{len(aa_seq) * 3} != {len(result)}')
    # # fill with trailing gaps so each induced dna sequence is of the same length
    # result += (len(aa_seq)*3-len(result))*'-'
    return result


def induce_msa(logger, aa_msa_path, dna_ms_path, output_path):
    gene_name_to_aligned_aa_sequence = load_header2sequences_dict(aa_msa_path)
    gene_name_to_unaligned_dna_sequence = load_header2sequences_dict(dna_ms_path)

    result = ''
    with open(dna_ms_path) as f:
        for line in f:
            if line.startswith('>'):
                gene_name = line.lstrip('>').rstrip()
                induced_dna_sequence = induce_sequence(logger, gene_name_to_aligned_aa_sequence[gene_name], gene_name_to_unaligned_dna_sequence[gene_name])
                result += f'>{gene_name}\n{induced_dna_sequence}\n'

    with open(output_path, 'w') as f:
        f.write(result)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('aa_msa_path', help='path to a file with aligned proteins')
    parser.add_argument('dna_ms_path', help='path to a file with unaligned dna sequences')
    parser.add_argument('output_path', help='path to a file in which the induced dna alignment will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        induce_msa(logger, args.aa_msa_path, args.dna_ms_path, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
