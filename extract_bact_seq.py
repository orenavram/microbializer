# -*- coding: utf-8 -*-
"""
Created on Tue May 22 14:52:18 2018

@author: shirportugez

This script is a module that contains functions for extracting sequences.

"""
#Imports
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
#Functions
#******************************************************************************************************        
def read_sequence_file_to_dict(file_path, file_format = "fasta"):
    '''
    This function takes a path to sequence file as an input and reads the file into a dictionary using the Bio.SeqIO.to_dict()
    '''
    recs = SeqIO.parse(file_path, file_format)
    return SeqIO.to_dict(recs)    
#******************************************************************************************************        
def get_strand(start, end):
    '''
    This function takes start and end BLAST coordinates as input and returns the strand:
        '+' -> if start < end
        '-' -> if end >= start
    '''
    #Plus strand
    if int(start) < int(end):
        return '+'
    #Minus strand
    else:
        return '-'
#******************************************************************************************************            
def trim_SeqRecord_by_coords_and_strand(full_rec, start, end, strand):
    '''
    This function gets a SeqRecord and returns a new SeqRecord according to start and end coordinates and a strand.
    '''
    start, end = int(start), int(end)
    #Extract sequence, get reverse complement if the strand is '-'
    if strand == '+':
        rec = full_rec[start - 1 :end]     
    else:
        rec = full_rec[end - 1 :start].reverse_complement()
    
    return rec
    
#******************************************************************************************************            
def extract_seq(accession, start, end, seqs_dict):
    '''
    This function extracts sequence according to accesion number and blast start and end coordinates.
    The function extracts the sequence from sequences dict (that needs to be generated using the function read_sequence_file_to_dict). 
    '''
    
    #Get strand
    strand = get_strand(start, end)
    
    #Get the SeqRecord (full sequence) of the accession number from the input fasta file 
    rec = seqs_dict[accession]
    
    #trim sequence
    try:
        trimed_rec = trim_SeqRecord_by_coords_and_strand(rec, start, end, strand)
    except (TypeError):
        print("Problem with sequence extraction.\n Please make sure the specified accession number is correct and exists in the input sequence file!")
    
    trimed_rec.id = "\t".join([accession,str(start),str(end)])
    trimed_rec.description = ""
    
    return trimed_rec   
#******************************************************************************************************        
def read_bact_genes_file(genes_file_path, seqs_file, output_seqs_path, seqs_file_format = "fasta"):
    '''
    This function gets a path to specific bacterial genes file (see explanation below), sequences file 
    (should be genbank or fasta), sequences file format (genbank or fasta) and to an output file. 
    The function creates new fasta output file with all the gene sequences from the input 
    file (i.e. sequence for each raw from the input file). 
    
    The input genes file should be a tab delimited file of the following form:
        -------------------------------------------------------
        |gene1_accesion_number    gene1_start    gene1_end    |
        |gene2_accesion_number    gene2_start    gene2_end    |
        |gene3_accesion_number    gene3_start    gene3_end    |
        -------------------------------------------------------
        
    '''
    
    #Read sequence input file into a dictionary
    seqs_dict = read_sequence_file_to_dict(seqs_file, seqs_file_format)
    
    #Read input genes file
    df = pd.read_csv(genes_file_path, sep = "\t", header = None)
    assert(len(df.columns) == 3) #maske sure there are three columns in the genes input file
    df.columns = ["accesion_number", "start", "end"]
    
    #Extract sequence of each row
    df_recs = df.apply(lambda x: extract_seq(x["accesion_number"], x["start"], x["end"], seqs_dict), axis = 1)

    #Create output fasta file
    try:
        SeqIO.write(list(df_recs), output_seqs_path, "fasta") #Convert records series to list (for output purposes)
    except:
        print("Problem with output sequence file writing!")
    
#******************************************************************************************************          



if __name__ == '__main__':
    try:
        genes_file_path, seqs_file, output_seqs_path = sys.argv[1:4]
    except (ValueError,IndexError):
        print("The following arguments should be given:\n1)path to genes file\n2)path to sequence file (fasta or genbank)\n3)path to output file\n")
    
    if len(sys.argv) == 5:
        seqs_file_format = sys.argv[4]
        read_bact_genes_file(genes_file_path, seqs_file, output_seqs_path, seqs_file_format)
    else:
        read_bact_genes_file(genes_file_path, seqs_file, output_seqs_path)
    