# -*- coding: utf-8 -*-
"""
Created on Wed May 23 15:53:14 2018

@author: shirportugez

This module creates input files for mcl & runs mcl.

Input:
    **blast dir - dir to blast result files
    **OGs csv file- two dimentional tab delimited file with OGs as rows and bacteria as columns. 
    Each cell will contain the accession_start_end coordinates of the gene (we might want to change it to geneID in the future, but then, I'm not sure how we will represent unannotated genes).
    
    example:
       --------------------------
       |   |bact1 |bact2 |bact3 |
       |OG1|ac-s-e|ac-s-e|ac-s-e|
       |OG2|ac-s-e|ac-s-e|ac-s-e|
       |OG3|ac-s-e|ac-s-e|ac-s-e|
       --------------------------
    **mcl dir - path to existing/unexisting folder that will contain the output orthologs table, 
    new subdirectories named input_mcl and output_mcl will be created (will contain the input files
    and output files for mcl respectively).
Output:
    orthologs csv file named ortholog_groups_table.csv, under the mcl_dir.
    This file will contain all ortholog groups from the og input file 
    (i.e. only OGs that were clustered into a single cluster).
"""

#Imports
import sys
import os
import pandas as pd
from itertools import combinations
#Functions
#******************************************************************************************************  
def get_all_pair_genes_combinations_of_specific_OG(row):
    '''
    This function gets a single row from the OGs file and returns a list of lists 
    (each inner list contains the og, followed by the first and second genes (the combination)). 
    '''
    og = row['OG'] #Get the OG label
    OGs_genes = row.drop(labels = "OG", axis = 0).values #Get all genes 
    combos_l = list(combinations(OGs_genes, 2)) #Get combinations of all pairs of genes
    combos_l = [list(t) for t in combos_l]
    [e.insert(0,og) for e in combos_l]
    return combos_l
#******************************************************************************************************  
def list_of_lists_to_flaten_list(ll):
    '''
    This function takes as input a list of lists and returns a flaten list 
    (single list).
    Example: [[1,2],[3,4]] --> [1,2,3,4]
    '''
    return [e for l in ll for e in l] #Extract series values into list of lists            
#******************************************************************************************************  
def get_genes_accs_start_end(gene, sep = "-"):
    '''
    This function splits a string representing a gene: accesion, start, end that are seperated by specific delimiter (default is "-")
    '''
    return gene.split(sep)
#******************************************************************************************************  
def get_blastres_path_by_two_acces(blast_d, acc1, acc2):
    '''
    This function gets a path to blast results folder that contains multiple blast result files and extracts the blast results path of
    the blast of acc1 vs acc2 (acc1 as a query, acc2 as a subject)
    '''
    files = os.listdir(blast_d) #List of all files under the blast dir
    query = "_".join([acc1,acc2]) #File name to search for - gene1 as query, gene2 as subject
    file = [f for f in files if query in f] #Get file names that match 
    assert(len(file) == 1)
    file_path = os.path.join(blast_d, file[0])
    return file_path
#******************************************************************************************************  
def extract_similarity_btw_two_genes(x, blast_d, score = "bitscore"):
    '''
    This function extracts the similarity between two genes (1,2) from a blast results file. The similarity is defined as
    their bit score (at least for now =] ).
    Input: path to blast dir with all blast result files (we will look for a file named bact1_bact2 in this dir, when bact1-query, bact2-subject), "gene1" (ac_s_e of gene1) "gene2" (ac_s_e of gene 2), score (set to bitscore for as default).
    Output: number representing the similarity of gene 1 and gene 2 (bitscore for now).
    '''
    g1, g2 = x['gene1'], x['gene2'] #REMOVE- toy example
    
    #Get gene1, gene2 accesion, start and end coordinates
    acc1, s1, e1 = get_genes_accs_start_end(g1)
    acc2, s2, e2  = get_genes_accs_start_end(g2)


    #Find blast file and read blast file (headers are included)
    blast_p = get_blastres_path_by_two_acces(blast_d, acc1, acc2)
    blast_df = pd.read_csv(blast_p, sep = "\t", header = 0)   
    
    #Get score from blast results - get column with g1 and g2
    mask = (blast_df['qseqid'] == acc1) & (blast_df['qstart'].astype(str) == s1) & (blast_df['qend'].astype(str) == e1) & (blast_df['sseqid'] == acc2) & (blast_df['sstart'].astype(str).isin([s2,e2])) & (blast_df['send'].astype(str).isin([s2,e2]))
    r = blast_df[mask] #Get the row
    assert(len(r) == 1) #Make sure the returned df has only a single row
    return r[score].values[0] #Return score
#******************************************************************************************************      
def create_mcl_input_file(x, mcl_input_d):
    '''
    This function gets a df with raws representing all gene combination of specific OG and creates an input
    file for mcl for this OG.
    '''
    og = str(x['OG'].values[0]) #Get og
    og_mcl_p = os.path.join(mcl_input_d, "og" + og + ".mcl_input") #path to mcl file of this og
    out_df = x.drop('OG', axis = 1)
    out_df.to_csv(og_mcl_p, sep = "\t", header = False, index = False)
    return og_mcl_p
#******************************************************************************************************      
def run_mcl(x, threshold = None):
    '''
    This function runs mcl, using an input mcl file and an output. Threshold variable is optional
    in case we want to exclude edges below specific weight.
    '''
    input_p, output_p = x["input_path"], x["output_path"]
    if threshold: #Threshold is specified
        cmd = "mcl {} --abc -tf 'gq({})' -o {}".format(input_p, threshold, output_p)
    else: #No threshold specified
        cmd = "mcl {} --abc -o {}".format(input_p, output_p)
    os.system(cmd)
#******************************************************************************************************      
def bool_is_file_contain_i_lines(x, i):
    '''
    This function gets a path to a file and an integer and returns True if the 
    file contains i lines, False otherwise.
    '''
    with open(x,"r") as f:
        lines = f.read()
    if lines.count("\n") == i:
        return True
    else:
        return False
#******************************************************************************************************      
if __name__ == '__main__':
    try:
        blast_dir, og_path, mcl_dir = sys.argv[1:4]
    except (ValueError,IndexError):
        print("The following arguments should be given:\n1)dir containing blast results files\n2)OG file (bacteria as columns, OGs as rows)\n3)MCL dir for input and output files\n")

    #Create mcl input & output dirs to hold the input and output mcl files
    mcl_input_d, mcl_output_d = os.path.join(mcl_dir,"input_mcl") , os.path.join(mcl_dir,"output_mcl") 
    os.makedirs(mcl_dir, exist_ok=True) #Create dir
    os.makedirs(mcl_input_d, exist_ok=True) #Create dir
    os.makedirs(mcl_output_d, exist_ok=True) #Create dir

    #Read input OGs file (see description at the beginning of the page)
    df = pd.read_csv(og_path)
    df.columns = ['OG'] + list(df.columns[1:])
    
    #Get all genes combinations for each row (i.e. all ac_s_e of each row) - create a new df with row for each combination
    #Generate a new df with 4 columns that will contain a row for each gene combination - "OG" (the OG group name of this combination) "gene1" (ac_s_e of gene1) "gene2" (ac_s_e of gene 2), "score" (the blast bitscore of gene1 against gene2)
    combos_s = df.apply(get_all_pair_genes_combinations_of_specific_OG, axis = 1)
    combos_l = list_of_lists_to_flaten_list(combos_s.values) #Extract series values into list of lists 
    combos_df = pd.DataFrame(combos_l) #Create df - new row for each genes combination
    combos_df.columns = ['OG','gene1','gene2']
    
    #Get the score of each combination
    combos_df['score'] = combos_df.apply(extract_similarity_btw_two_genes, args = (blast_dir,), axis = 1)
    
    #Each OG (groups of rows with same "OG" label will be an intput to the mcl algorithm)
    mcl_ps = combos_df.groupby('OG').apply(create_mcl_input_file, mcl_input_d)
    mcl_ps = mcl_ps.to_frame().rename(columns = {0:"input_path"}) 
    mcl_ps["output_path"] = mcl_ps["input_path"].apply(lambda x: x.replace("input","output"))
    
    #Run mcl on each input file
    mcl_ps.apply(run_mcl, axis = 1)
    
    #Find out which groups remain contact
    mcl_ps['single_cluster_flag'] = mcl_ps["output_path"].apply(bool_is_file_contain_i_lines, args = (1,))
    
    #Create final orthologs table 
    full_df = pd.concat([df, mcl_ps['single_cluster_flag']], axis = 1)
    df = full_df[full_df['single_cluster_flag']] 
    df = df.drop('single_cluster_flag', axis = 1)
    ortho_grps_p = os.path.join(mcl_dir,"ortholog_groups_table.csv")
    df.to_csv(ortho_grps_p, sep = ",", header = 0, index = False)









