# -*- coding: utf-8 -*-
"""
Created on Tue May 29 20:27:58 2018

@author: student
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 22 14:52:18 2018

@author: shirportugez

This script is a module that extract sequences of ortholog groups. 
Input: 
    (1) a path to ortholkogs table.
    (2) directory for sequence files (each bacteria in the og table should have one sequence file with the name of the bacteria that should contain all genes of that bacteria).
    (3) a directory for output sequence files - if the folder doesn't exists, the script will create it.
    (4) input and output sequence file foramt (fasta or genebank)
Output:
    directory containing fasta files for each ortholog group.
"""

#Imports
import sys
import os
from src.constract_final_table import list_of_lists_to_flaten_list
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
def get_og_genes(row):
    '''
    This function gets a single row from the OGs file and returns a list of lists: each inner list contain the og as the 
    first element and the gene name as the second.
    '''
    og = row['OG'] #Get the OG label
    OGs_genes = row.drop(labels = "OG", axis = 0) #Get all genes 
    genes_l = []
    [genes_l.append([og,bact,gene]) for bact,gene in zip(OGs_genes.index , OGs_genes.values)] #Create list of lists ([['og',bacteriaa, gene_name1],[og,bacteria2, gene_name2],...])
    return genes_l
#******************************************************************************************************          
def get_bact_sequence_path(x, files_dict):
    '''
    This function gets bacteria name and files dictionary (file base name as a key, full path as a value) and returns file name that contains the bacterial name.
    '''
    bact = x['bacteria']
    bact_files = [f for f in files_dict.keys() if bact == f]
    assert(len(bact_files) ==1) #Only single file for this bacteria
    return files_dict[bact_files[0]]
#******************************************************************************************************          
def extract_gene_sequence(x, file_format):
    '''
    This function gets a row from the og_genes_and_paths_df and file format ("fasta" or "genbank" for SeqIO package), and returns a string representing the sequence of the 
    gene.
    '''
    recs = read_sequence_file_to_dict(x['seq_path'],file_format)
    assert(x['gene'] in recs.keys())
    return recs[x['gene']]    
#******************************************************************************************************              
def write_seqfile(x, file_format):
    '''
    This function gets a row from og_recs ( that contains 'SeqRecords' column that contains list of 
    Seqrecords, 'fasta_output_path' that contains the path to the output file and file format) and 
    creates new sequence file in the required format.
    '''
    SeqIO.write(x['SeqRecords'], x['fasta_output_path'], file_format)
#******************************************************************************************************              
if __name__ == '__main__':
    try:
        og_path, seqs_dir, output_seqs_d, seqs_file_format = sys.argv[1:]
    except (ValueError,IndexError):
        print("The following arguments should be given:\n1) a path to final orthologs table\n2)a path to directory that contains all bacteria sequence files (fasta file for each bacteria species)\n3)path to output directory for fasta files of each OG\n")
    
    #Create output directory if doesn't exists - will include all OGs fasta files
    os.makedirs(output_seqs_d, exist_ok=True) #Create dir
    
    #Read input OGs file (see description at the beginning of the page)
    df = pd.read_csv(og_path, header = 0)
    
    #Get df with the path to each bacteria species 
    files = os.listdir(seqs_dir)
    files_base = [os.path.splitext(f)[0] for f in files] #Files without extensions
    full_files = [os.path.join(seqs_dir, f) for f in files]
    files_dict = dict(zip(files_base, full_files))
    bact_seqs_df = pd.DataFrame(data = list(df.columns[1:]), columns = ['bacteria'], index = range(len(df.columns[1:]))) #df with bacteria species as index, and column of sequence_path contains a string with the path to the specific bacteria sequence file
    bact_seqs_df['seq_path'] = bact_seqs_df.apply(get_bact_sequence_path, args = (files_dict,), axis = 1)
    
    #Get all genes for each row of the df
    og_genes_s = df.apply(get_og_genes, axis = 1) #Get all genes, each series value is a list of lists [['og',bact1, gene1],[og,bact2,gene2],...]
    og_genes_l = list_of_lists_to_flaten_list(og_genes_s.values)
    og_genes_df = pd.DataFrame(og_genes_l)
    og_genes_df.columns = ['OG','bacteria','gene']
    
    #Extract the sequence of each gene from the og genes df
    og_genes_and_paths_df = pd.merge(og_genes_df,bact_seqs_df, on = 'bacteria')
    og_genes_df['sequence'] = og_genes_and_paths_df.apply(extract_gene_sequence, args = (seqs_file_format,), axis = 1)

    #Create fasta file for each ortholog group
    og_recs = og_genes_df.groupby('OG').apply(lambda x: list(x['sequence'].values)).to_frame().reset_index().rename(columns = {0:'SeqRecords'}) #Get df with SeqRecords list for each OG
    og_recs['fasta_output_path'] = og_recs['OG'].apply(lambda x: os.path.join(output_seqs_d,"og_" + str(x) + ".fa")) #Get output path for each OG
    og_recs.apply(write_seqfile, args = (seqs_file_format,), axis = 1)
    