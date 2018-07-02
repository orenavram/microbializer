import subprocess
import pandas as pd

def create_blast_DB(reference_seq, dbtype, output_path):
    '''
    input:  sequnce to base the DB on
            DB type (nucl/prot)
            path to output file
    output: balst DB based on the reference sequence
    '''
    cmd = 'makeblastdb -in {seq} -dbtype {db_type} -out {reference_DB}'.format(seq = reference_seq, db_type = dbtype, reference_DB = output_path)
    subprocess.call(cmd.split(), shell = True)
        

def blast_all_vs_all(program, query_seq, reference_DB, output_path):
    '''
    input:  type of blast to run(blastn/blastp)
            fasta file of query genome
            blast DB file
    output: query_vs_reference blast results file
    '''
    cmd = ('{prog} -query {query} -db {reference} -out {output} -max_target_seqs 2 -outfmt 6'.format(prog = program, query = query_seq, reference = reference_DB, output = output_path))
    subprocess.call(cmd.split(), shell = True)

def filter_blast(query_vs_reference, precent_identity_cutoff, e_value_cutoff, bit_score_cutoff, output_path):
    '''
    input:  path to file of blast results
            desiered cutoff values
    output: file with filtered results
    '''
    header = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    df = pd.read_csv(query_vs_reference, header = None, sep = '\t', names = header.split())
    result = df[(df.pident>=precent_identity_cutoff)  & (df.evalue<=e_value_cutoff) & (df.bitscore>=bit_score_cutoff)]
    result.to_csv(output_path, sep = '\t', index = False)  

if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    if len(argv) < 5:
        logger.error('Usage: blast_all_vs_all <program> <query_sequence_file> <reference_seq_file> <output_path>')
        exit()
    else:
        blast_all_vs_all(*argv[1:])
