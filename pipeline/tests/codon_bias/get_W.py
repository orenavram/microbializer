import subprocess
import sys
from sys import argv
import argparse
import logging
import os
from Bio import SeqIO
import CodonUsageModified as CodonUsage
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries import consts

HITS_TO_KEEP_FOR_EACH_REFERENCE_HEG = 3


def create_blast_db(ORF_file_path, logger):
    """
    input: file path to nucleotide fasta of open reading frames

    output: nucleotide blast DB based
    """
    output_file = ORF_file_path + '.db'
    dbtype = "nucl"
    cmd = 'makeblastdb -in {} -out {} -dbtype {}'.format(ORF_file_path, output_file, dbtype)
    logger.info('Making blastdb with command: ' + cmd)
    subprocess.run(cmd, shell=True)


def blast_with_HEG(ORF_file_path, logger):
    """
    input: file path to nucleotide fasta of open reading frames

    output: 40 highest match hits from a tblastn with the ecoli reference query
    """
    output_file = ORF_file_path + "_HEG_hits"
    database = ORF_file_path + '.db'

    cmd = f'tblastn -db {database} -query {consts.HEGS_ECOLI_FILE_PATH} -out {output_file} -outfmt 6 -max_target_seqs {HITS_TO_KEEP_FOR_EACH_REFERENCE_HEG}'
    logger.info("Finding Hits with command: \n" + cmd)
    subprocess.run(cmd, shell=True)

    # Delete Blastdb
    cmd = f'rm {ORF_file_path}.*'
    logger.info("Cleaning directory with command: " + cmd)
    subprocess.run(cmd, shell=True)


def get_hits_only(ORF_file_path, logger):
    """
    input: file path to nucleotide fasta of open reading frames

    output: isolate the sequence ids of the hits
    """
    # Get sequence names in a temporary file
    blastfile = ORF_file_path + "_HEG_hits"
    hits_only_file = ORF_file_path + "_HEG_hits_only.txt"
    cmd = f'cut -f 2 {blastfile} | sort -u > {hits_only_file}'
    subprocess.run(cmd, shell = True)
    logger.info("getting sequence names into "+ hits_only_file)

    # Remove Blast Output File
    cmd = f'rm {blastfile}'
    subprocess.run(cmd, shell = True)
    logger.info("removing blast file " + blastfile)


def make_HEF_fa(ORF_file_path, logger, output_file):
    """
    input: file path to nucleotide fasta of open reading frames

    output: fasta file containg the selected highly expressed genes
    """
    #Create output file
    temp_output_file = open(ORF_file_path + "_HEG.fa", "w")
    hits_only_file = ORF_file_path + "_HEG_hits_only.txt"
    logger.info("Creating Outputfile ")

    #Parse Sequences with
    tempfile = open(hits_only_file, "r")
    hits = tempfile.readlines()
    for hit in hits:
        with open(ORF_file_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if(hit[0:-1] == record.id):
                   # if(record.seq.find('N'))
                    temp_output_file.write(record.format("fasta"))
                    break

    #Move hits only file
    if not os.path.exists(output_file):
        os.makedirs(output_file+"/HEG_hits")
    cmd = f'mv {hits_only_file} {output_file}/HEG_hits'
    subprocess.run(cmd, shell = True)
    logger.info("moving temporary file "+ hits_only_file)


def save_index(output_file, genomeIndex, filename, filepath):

    filepath = (output_file + "/genomeIndex")

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    file = filepath + "/" + filename[0:-8]

    with open(file, "w") as file:
        json.dump(genomeIndex.index, file)


def calculate_codon_bias(filepath, filename, logger, output_file):
    """
    input: file path to nucleotide fasta of the highly expressed genes
           name of the ORF file

    output: fasta file containg the selected highly expressed genes
    """

    genomeIndex = CodonUsage.CodonAdaptationIndex()
    genomeIndex.generate_index(filepath + "_HEG.fa")

    save_index(output_file, genomeIndex, filename, filepath)

    cmd = f'rm {filepath}_HEG.fa'
    subprocess.run(cmd, shell = True)
    logger.info("removing HEF fasta file "+ filepath + "_HEG.fa")


def get_W(ORFs_file, output_file, logger):
    logger.info("Getting HEGs from the ORFs file: " + str(ORFs_file))

    # Create Blast Database from the ORFs file
    create_blast_db(ORFs_file, logger)

    # Identify HEGs in the ORFs file
    blast_with_HEG(ORFs_file, logger)
    get_hits_only(ORFs_file, logger)
    make_HEF_fa(ORFs_file, logger, output_file)

    # Find Codon bias of HEGs
    calculate_codon_bias(filepath, filename, logger, output_file)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('ORF_file', help='path to input ORF file')
    parser.add_argument('output_dir', help='path to output dir')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        get_W(args.ORF_file, args.output_dir, logger)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
