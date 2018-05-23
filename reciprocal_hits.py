import regex as re
import pandas as pd
import csv


def reciprocalHits(blast_out1, blast_out2, output_path):
    """
    :param blast_out1: BLAST output file comparing org1 with org2 delimitered by any space charecter
    :param blast_out2: BLAST output file comparing org2 with org1 delimitered by any space charecter
    creating csv file containng reciprocal hits of balst_out1 & blast_out2
    """
    gene_pairs = generatePairingDict(blast_out1)  #creating dictionary of genes pairs of org_1 compared with org_2
    index = 0
    lines_counter = 0
    source_file2 = open(blast_out2, "r")
    lines = source_file2.readlines()
    for line in lines:
        if lines_counter == 0: #getting organisms name from first line in the BLAST file
            org_1, org_2 = getOrganismsName(line)
            df = pd.DataFrame(columns=[org_1, org_2]) #creating data frame containing coulmn for each organism.
        else: #checking if genes pairs of org_2 compared with org_1 are found in org_1 compared with org_2
            gene, ortho_gene = getGenesPair(line)
            if ortho_gene in gene_pairs:
                if gene == gene_pairs[ortho_gene ]:
                    df.loc[index]=[gene,ortho_gene ]
                    index += 1
        lines_counter += 1
    source_file2.close()
    df.to_csv(output_path + "/reciprocal_hits_ " + org_1+" _" +org_2 + ".csv", index=False) #pushing data frame to csv file


def getOrganismsName(blast_out_line):
        """
        :param blast_out_line: first line of BLAST output file comparing org1 with org2
        :return (organism_1, organism_2): names of two organisms compared
        """
        res = re.search("(\w+)(\s{1})(\w+)",blast_out_line)  # look for string containing 2 genes separated by a space character
        if res:
            org_1 = re.search("(\w+)(\s{1})(\w+)", blast_out_line).group(1);
            org_2 = re.search("(\w+)(\s{1})(\w+)", blast_out_line).group(3);
        return (org_1, org_2)


def getGenesPair(blast_out_line):
    """
    :param blast_out_line: line of BLAST output file comparing org1 with org2
    :return (gene, ortho_gene): pair of orthologous genes.
    """
    res = re.search("(\w+)(\s{1})(\w+)", blast_out_line)
    if res:
        gene = re.search("(\w+)(\s{1})(\w+)", blast_out_line).group(1);
        ortho_gene = re.search("(\w+)(\s{1})(\w+)", blast_out_line).group(3);
        return (gene, ortho_gene)

def generatePairingDict(blast_out):
    """
    :param blast_out: BLAST output file comparing org1 with org2
    :return gene_pairs_dict: dictionary of genes as keys and matching ortho genes as values
    """
    gene_pairs_dict = {};
    source_file = open(blast_out, "r")
    lines = source_file.readlines()
    for line in lines:
        gene = re.search("(\w+)(\s{1})(\w+)", line).group(1);
        ortho_gene = re.search("(\w+)(\s{1})(\w+)", line).group(3);
        gene_pairs_dict[gene] = ortho_gene;
    source_file.close()
    return gene_pairs_dict;