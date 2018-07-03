import subprocess
import pandas as pd


def create_blast_DB(reference_seq, dbtype, output_path):
    '''
    input:  sequnce to base the DB on
            DB type (nucl/prot)
            path to output file
    output: balst DB based on the reference sequence
    '''
    cmd = 'makeblastdb -in {seq} -dbtype {db_type} -out {reference_DB}'.format(seq=reference_seq, db_type=dbtype,
                                                                               reference_DB=output_path)
    subprocess.call(cmd.split(), shell=True)


def blast_all_vs_all(program, query_seq, reference_DB, output_path):
    '''
    input:  type of blast to run(blastn/blastp)
            fasta file of query genome
            blast DB file
    output: query_vs_reference blast results file
    '''
    cmd = ('{prog} -query {query} -db {reference} -out {output} -outfmt 6'.format(prog=program, query=query_seq,
                                                                                  reference=reference_DB,
                                                                                  output=output_path))
    subprocess.call(cmd.split(), shell=True)
