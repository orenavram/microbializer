import sys
import os
from sys import argv
import argparse
import logging
import traceback
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


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
    gene_name_to_aligned_aa_sequence = {record.id: record.seq for record in SeqIO.parse(aa_msa_path, 'fasta')}
    gene_name_to_unaligned_dna_sequence = {record.id: record.seq for record in SeqIO.parse(dna_ms_path, 'fasta')}

    induced_dna_records = []
    for gene_name, dna_sequence in gene_name_to_unaligned_dna_sequence.items():
        induced_dna_sequence = induce_sequence(logger, gene_name_to_aligned_aa_sequence[gene_name], dna_sequence)
        induced_dna_records.append(SeqRecord(id=gene_name, name=gene_name, seq=Seq(induced_dna_sequence)))

    SeqIO.write(induced_dna_records, output_path, 'fasta')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('aa_msa_path', help='path to a file with aligned proteins')
    parser.add_argument('dna_ms_path', help='path to a file with unaligned dna sequences')
    parser.add_argument('output_path', help='path to a file in which the induced dna alignment will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        induce_msa(logger, args.aa_msa_path, args.dna_ms_path, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
