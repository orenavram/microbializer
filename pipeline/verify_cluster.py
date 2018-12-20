# -*- coding: utf-8 -*-
"""
Running example:
script_name.py /Users/Oren/Dropbox/Projects/microbializer/mock_output/all_recip_hits.csv /Users/Oren/Dropbox/Projects/microbializer/mock_output/splitted_putative_orthologs_table/NC_003197.2_cds_NP_459006.1_1_6.splitted_putative_orthologs_table /Users/Oren/Dropbox/Projects/microbializer/mock_output/tmp_dir/verified_clusters /Users/Oren/Dropbox/Projects/microbializer/mock_output/verified_clusters/NC_003197.2_cds_NP_459006.1_1_6.verified_clusters

Created on Wed May 30 14:29:48 2018

@author: shirportugez

This module gets an putative orthologs table, runs mcl on each ortholog group and returns new final orthologs table
with OGs that .

Input:
    ** path to reciptocal hits file
        ------------------------
       |g1\tg2\tsimilarity score|
       |g2\tg3\tsimilarity score|
       |g1\tg4\tsimilarity score|
       -------------------------

    **path to OGs csv file- two dimentional tab delimited file with OGs as rows and bacteria as columns.

       --------------------------
       |   |bact1 |bact2 |bact3 |
       |OG1|ac-s-e|ac-s-e|ac-s-e|
       |OG2|ac-s-e|ac-s-e|ac-s-e|
       |OG3|ac-s-e|ac-s-e|ac-s-e|
       --------------------------
    **tmp dir - path to tmp folder that will contain input and output files for mcl. 
    They both will be removed after the function is done.
    ** output path 
    
Output:
    orthologs csv file with rows from the input OG table that are orthologs (all genes have been clustered together).
    
"""

import sys
import os
import pandas as pd
from itertools import combinations

#******************************************************************************************************
def get_all_pair_genes_combinations_of_specific_OG(row):
    '''
    This function gets a single row from the OGs file and returns a list of lists 
    (each inner list contains the og, followed by the first and second genes (the combination)). 
    '''
    og = row['OG_name'] #Get the OG label
    OGs_genes = row.drop(labels = 'OG_name', axis = 0).values #Get all genes
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
def extract_similarity_btw_two_genes(x, reciphits_df):
    '''
    This function extracts the similarity between two genes (1,2) from a reciprocal hits file. 
    Input: row from gene combination table (df row that contains g1 and g2), df of to reciprocal hits 
    Output: number representing the similarity of gene 1 and gene 2 (bitscore for now).
    '''
    g1, g2 = x['gene1'], x['gene2'] 
    
    #Get score from blast results - get column with g1 and g2
    mask = ((reciphits_df['gene1'] == g1) & (reciphits_df['gene2'] == g2)) | ((reciphits_df['gene2'] == g1) & (reciphits_df['gene1'] == g2))
    r = reciphits_df[mask] #Get the row
    assert(len(r) == 1) #Make sure the returned df has only a single row
    return r['similarity'].values[0] #Return score
#******************************************************************************************************      
def create_mcl_input_file(x, mcl_input_d):
    '''
    This function gets a df with rows representing all gene combination of specific OG and creates an input
    file for mcl for this OG.
    '''
    og = x['OG_name'].values[0] #Get og name
    og_mcl_p = os.path.join(mcl_input_d, og + ".mcl_input") #path to mcl file of this og
    out_df = x.drop('OG_name', axis = 1)
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
def bool_is_file_contain_i_lines(x):
    '''
    This function gets a path to a file and an integer and returns True iff the
    file contains a single line (to check how many clusters mcl found).
    '''
    with open(x) as f:
        i=0
        for line in f:
            i+=1
            if i>1:
                return False
        return i == 1

#******************************************************************************************************      
if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    logger.info(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    recip_hits_path, og_path, tmp_dir, out_og_path = sys.argv[1:]
    #except (ValueError,IndexError):
    #    logger.error("The following arguments should be given:\n1)reciprocal hits file\n2)putative OG file (bacteria as columns, OGs as rows)\n3)tmp dir for mcl input and output files\n4)output path for final orhoglogs table ")

    #Create tmp dir if doesn't exists - and subdirectories fro mcl input and output files
    os.makedirs(tmp_dir, exist_ok=True) #Create dir
    mcl_input_d, mcl_output_d = os.path.join(tmp_dir,"input_mcl") , os.path.join(tmp_dir,"output_mcl") 
    os.makedirs(mcl_input_d, exist_ok=True) #Create dir
    os.makedirs(mcl_output_d, exist_ok=True) #Create dir

    #Read input OGs file (see description at the beginning of the page)
    df = pd.read_csv(og_path, header=0) #TODO: make sure this change is OK
    df.rename(columns={df.columns[0]: 'OG_name'}, inplace=True) #TODO: remove GENE label and use OG

    #Read reciprocal blast hits table
    reciphits_df = pd.read_csv(recip_hits_path, header = None)
    reciphits_df.columns = ['gene1', 'gene2', 'similarity']
        
    #Get all genes combinations for each row (i.e. all ac_s_e of each row) - create a new df with row for each combination
    #Generate a new df with 4 columns that will contain a row for each gene combination - "OG" (the OG group name of this combination) "gene1" (ac_s_e of gene1) "gene2" (ac_s_e of gene 2), "score" (the blast bitscore of gene1 against gene2)
    combos_s = df.apply(get_all_pair_genes_combinations_of_specific_OG, axis = 1)
    combos_l = list_of_lists_to_flaten_list(combos_s.values) #Extract series values into list of lists 
    combos_df = pd.DataFrame(combos_l) #Create df - new row for each genes combination
    combos_df.columns = ['OG_name','gene1','gene2']
    combos_df.dropna(inplace = True)
	
    #Get the score of each combination
    combos_df['score'] = combos_df.apply(extract_similarity_btw_two_genes, args=(reciphits_df,), axis = 1)
    
    #Each OG (groups of rows with same "OG" label will be an input to the mcl algorithm)
    mcl_ps_series = combos_df.groupby('OG_name').apply(create_mcl_input_file, mcl_input_d)
    mcl_ps = mcl_ps_series.to_frame().rename(columns = {0:"input_path"})
    mcl_ps['output_path'] = mcl_ps['input_path'].apply(lambda x: x.replace('input', 'output'))
    
    #Run mcl on each input file
    mcl_ps.apply(run_mcl, axis=1)

    #Find out which groups remain contact
    mcl_ps['single_cluster_flag'] = mcl_ps['output_path'].apply(bool_is_file_contain_i_lines)
    
    #Create final orthologs table 
    full_df = pd.concat([df.set_index('OG_name'), mcl_ps['single_cluster_flag']], axis=1)
    df = full_df[full_df['single_cluster_flag']] 
    df.drop('single_cluster_flag', axis=1, inplace=True)
    df.to_csv(out_og_path, sep=',', header=0, index=False)









