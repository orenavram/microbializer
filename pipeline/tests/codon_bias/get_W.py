import subprocess
import sys
from sys import argv
import argparse
import logging
import os
from Bio import SeqIO
import CodonUsageModified as CodonUsage
import pipeline_auxiliaries as pa
import json

def create_blast_DB(ORF_file_path, logger):
    """
    input: file path to nucleotide fasta of open reading frames

    output: nucleotide blast DB based on the reference sequence
    """
    output_file = ORF_file_path + '.db'
    dbtype = "nucl"
    cmd = 'makeblastdb -in {} -out {} -dbtype {}'.format(ORF_file_path, output_file, dbtype)
    logger.info('Making blastdb with command: ' + cmd)
    subprocess.run(cmd, shell = True)


def blast_with_HEG(ORF_file_path, HEG_reference_file, logger):
    """
    input: file path to nucleotide fasta of open reading frames

    output: 40 highest match hits from a tblastn with the ecoli reference query
    """
    output_file = ORF_file_path + "_HEG_hits"
    database = ORF_file_path + '.db'
    NUMBER_BLAST_HITS = 3

    cmd = f'tblastn -db {database} -query {HEG_reference_file} -out {output_file} -outfmt 6 -max_target_seqs {NUMBER_BLAST_HITS}'
    logger.info("Finding Hits with command: \n" + cmd)
    subprocess.run(cmd, shell = True)

    #Delete Blastdb
    cmd = f'rm {ORF_file_path}.*'
    logger.info("Cleaning directory with command: " + cmd)
    subprocess.run(cmd, shell = True)



def get_hits_only(ORF_file_path, logger):
    """
    input: file path to nucleotide fasta of open reading frames

    output: isolate the sequence ids of the hits
    """
    #Get sequence names in a temporary file
    blastfile = ORF_file_path + "_HEG_hits"
    hits_only_file = ORF_file_path + "_HEG_hits_only.txt"
    cmd = f'cut -f 2 {blastfile} | sort -u > {hits_only_file}'
    subprocess.run(cmd, shell = True)
    logger.info("getting sequence names into "+ hits_only_file)

    #Remove Blast Output File
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

def get_W(ORF_dir, OG_dir, HEG_reference_file, output_file, logger, filename):
    logger.info("Getting HEGs from genome: "+ str(filename))

    filepath = ORF_dir + '/' + filename

    #Create Blast Database for each Genome
    create_blast_DB(filepath, logger)

    #Identify HEG in each genome
    blast_with_HEG(filepath, HEG_reference_file, logger)
    get_hits_only(filepath, logger)
    make_HEF_fa(filepath, logger, output_file)


    #Find Codon bias of HEGs
    calculate_codon_bias(filepath, filename, logger, output_file)

if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('ORF_dir', help='path to input fasta directory')
    parser.add_argument('OG_dir', help='path to input Orthologous group directory')
    parser.add_argument('HEG_reference_file', help='path to file of highly expressed bacterial genes')
    parser.add_argument('output_file', help= 'path to output location')
    parser.add_argument('file', help= 'Number of genome that will be calculated')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')


    args = parser.parse_args()

    level = logging.INFO
    logger = pa.get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        get_W(args.ORF_dir, args.OG_dir, args.HEG_reference_file, args.output_file, logger, args.file)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
