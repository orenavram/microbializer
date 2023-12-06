import subprocess
import sys
from sys import argv
import argparse
import logging
import os
import shutil
from Bio import SeqIO
import CodonUsageModified as CodonUsage
import json
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts

MAX_HITS_TO_KEEP_FOR_EACH_REFERENCE_HEG = 3
BLAST_IDENTITY_PERCENT_THRESHOLD = 35
BLAST_EVALUE_THRESHOLD = 0.01


def find_HEGs_in_orf_file(ORFs_file, genome_name, tmp_dir, logger):
    """
    input: file path to nucleotide fasta of open reading frames

    output: nucleotide blast DB based
    """

    # Create blast db from ORF file
    db_dir = os.path.join(tmp_dir, f'{genome_name}_db')
    db_name = os.path.join(db_dir, genome_name + '.db')
    cmd = f'makeblastdb -in {ORFs_file} -out {db_name} -dbtype nucl'
    logger.info('Making blastdb with command: \n' + cmd)
    subprocess.run(cmd, shell=True)

    # Query blast db with ecoli HEGs reference file
    hegs_hits_file = os.path.join(tmp_dir, genome_name + '_HEG_hits.tsv')
    cmd = f'tblastn -db {db_name} -query {consts.HEGS_ECOLI_FILE_PATH} -out {hegs_hits_file} -outfmt 6 ' \
          f'-max_target_seqs {MAX_HITS_TO_KEEP_FOR_EACH_REFERENCE_HEG}'
    logger.info("Finding Hits with command: \n" + cmd)
    subprocess.run(cmd, shell=True)

    shutil.rmtree(db_dir, ignore_errors=True)

    # Filter hits to find actual HEGs and write their names into a file
    hegs_df = pd.read_csv(hegs_hits_file, delimiter='\t', names=consts.BLAST_OUTPUT_HEADERS)
    hegs_df_filtered = hegs_df.loc[(hegs_df['identity_percent'] > BLAST_IDENTITY_PERCENT_THRESHOLD) & (hegs_df['evalue'] < BLAST_EVALUE_THRESHOLD)]
    hegs_names = set(hegs_df_filtered['subject'])
    HEGs_names_file_path = os.path.join(tmp_dir, genome_name + '_HEG_hits_only.txt')
    with open(HEGs_names_file_path, 'w') as HEGs_names_file:
        HEGs_names_file.write('\n'.join(hegs_names))
    logger.info("HEGs names were written into " + HEGs_names_file_path)

    return HEGs_names_file_path


def make_HEGs_fasta(ORFs_file, genome_name, HEGs_names_file_path, tmp_dir, logger):
    with open(HEGs_names_file_path, "r") as ORFs_hegs_file:
        HEGs_names = ORFs_hegs_file.read().splitlines()

    HEGs_fasta_file_content = ''
    with open(ORFs_file) as ORF_file:
        for record in SeqIO.parse(ORF_file, "fasta"):
            if record.id in HEGs_names:
                HEGs_fasta_file_content += record.format("fasta")

    HEGs_fasta_file_path = os.path.join(tmp_dir, genome_name + "_HEGs.fa")
    with open(HEGs_fasta_file_path, "w") as HEGs_fasta_file:
        HEGs_fasta_file.write(HEGs_fasta_file_content)

    logger.info(f'Created {HEGs_fasta_file_path} that contains the sequences of the HEGs')
    return HEGs_fasta_file_path


def calculate_codon_bias(HEGs_fasta_file_path, output_dir, genome_name, logger):
    genome_index = CodonUsage.CodonAdaptationIndex()
    genome_index.generate_index(HEGs_fasta_file_path)

    relative_adaptiveness_output_path = os.path.join(output_dir, genome_name + ".json")

    with open(relative_adaptiveness_output_path, "w") as relative_adaptiveness_file:
        json.dump(genome_index.index, relative_adaptiveness_file)

    logger.info(f'Calculated the relative adaptiveness of codons (W vector) for HEGs file {HEGs_fasta_file_path}, '
                f'and saved it in {relative_adaptiveness_output_path}')


def get_W(ORFs_file, output_dir, tmp_dir, logger):
    logger.info("Getting HEGs from the ORFs file: " + str(ORFs_file))

    genome_name = os.path.splitext(os.path.basename(ORFs_file))[0]

    # Identify HEGs in the ORFs file
    HEGs_names_file_path = find_HEGs_in_orf_file(ORFs_file, genome_name, tmp_dir, logger)
    HEGs_fasta_file_path = make_HEGs_fasta(ORFs_file, genome_name, HEGs_names_file_path, tmp_dir, logger)

    # Find Codon bias of HEGs
    calculate_codon_bias(HEGs_fasta_file_path, output_dir, genome_name, logger)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('ORF_file', help='path to input ORF file')
    parser.add_argument('output_dir', help='path to output dir')
    parser.add_argument('tmp_dir', help='path to tmp dir')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        get_W(args.ORF_file, args.output_dir, args.tmp_dir, logger)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
