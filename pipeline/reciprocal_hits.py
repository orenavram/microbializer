import regex as re
import pandas as pd
import os
import csv


def reciprocalHits(blast_out1, blast_out2, output_path, output_name):
    """
    :param blast_out1: BLAST output file comparing org1 with org2 delimitered by any space charecter
    :param blast_out2: BLAST output file comparing org2 with org1 delimitered by any space charecter
    creating csv file containng reciprocal hits of balst_out1 & blast_out2
    """
    gene_pairs = generatePairingDict(blast_out1)  #creating dictionary of genes pairs of org_1 compared with org_2
    index = 0
    lines_counter = 1
    source_file2 = open(blast_out2, "r")
    org_1 = (output_name.split('.')[0]).split('_vs_')[0]
    org_2 = (output_name.split('.')[0]).split('_vs_')[1]
    df = pd.DataFrame(columns=[org_1, org_2, 'bitscore'])
    for line in source_file2:
        gene, ortho_gene = getGenesPair(line)
        bitscore = getBitScore(line)
        if ortho_gene in gene_pairs:
            if gene == gene_pairs[ortho_gene]:
                df.loc[index]=[gene, ortho_gene, bitscore]
                index += 1
        lines_counter += 1
    source_file2.close()
    output_file = os.path.join(output_path, output_name)
    df.to_csv(output_file, index=False) #pushing data frame to csv file


def getOrganismsName(blast_out_line):
        """
        :param blast_out_line: first line of BLAST output file comparing org1 with org2
        :return (organism_1, organism_2): names of two organisms compared
        """
        res = re.split('\s+', blast_out_line)
        #org_1 = re.split('_cds_', res[0][4:-2])[0] #Shir - what will happend when we are not preddicting the cds for example? maybe the _cds assumption is not general enough?
        #org_2 = re.split('_cds_', res[1][4:-2])[0]
        org_1, org_2 = res[1], res[0]
        return (org_1, org_2)


def getGenesPair(blast_out_line):
    """
    :param blast_out_line: line of BLAST output file comparing org1 with org2
    :return (gene, ortho_gene): pair of orthologous genes.
    """
    res = re.split('\s+', blast_out_line)
    #gene = re.split('_cds_',res[0][4:-2])[1] #Shir - didn't understand this line, why to split according to _cds_? 
    #ortho_gene = re.split('_cds_',res[1][4:-2])[1] #Shir - didn't understand this line, why to split according to _cds_?
    gene, ortho_gene = res[1], res[0]
    return (gene, ortho_gene)

def getBitScore(blast_out_line):
    """
    :param blast_out_line: line of BLAST output file comparing org1 with org2
    :return bit_score of pair of genes compared in the line
    """
    res = re.split('\s+', blast_out_line.strip())
    bit_score = (res[-1])
    return bit_score

def generatePairingDict(blast_out):
    """
    :param blast_out: BLAST output file comparing org1 with org2
    :return gene_pairs_dict: dictionary of genes as keys and matching ortho genes as values
    """
    gene_pairs_dict = {}
    source_file = open(blast_out, "r")
    lines = source_file.readlines()
    lines_counter = 0
    for line in lines:
        if lines_counter > 0: #ignore first line
            gene, ortho_gene = getGenesPair(line)
            gene_pairs_dict[gene] = ortho_gene
        lines_counter += 1
    source_file.close()
    return gene_pairs_dict


if __name__ == '__main__':
        import logging
        import argparse
        logger = logging.getLogger('main')
        parser = argparse.ArgumentParser()

        parser.add_argument('blast_result1', help='path to blast result file of seq1 vs seq2')
        parser.add_argument('blast_result2', help='path to blast result file of seq2 vs seq1')
        parser.add_argument('output_path', help='path to output file')
        parser.add_argument('output_name', help='name of output file')
        args = parser.parse_args()

        reciprocalHits(args.blast_result1, args.blast_result2, args.output_path, args.output_name)
