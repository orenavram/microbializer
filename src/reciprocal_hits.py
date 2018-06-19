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
    for line in source_file2:
        if lines_counter == 1: #getting organisms name from first line in the BLAST file
            org_1, org_2 = getOrganismsName(line)
            df = pd.DataFrame(columns=[org_1, org_2, 'bitscore']) #creating data frame containing coulmn for each organism.
        elif lines_counter > 0: #checking if genes pairs of org_2 compared with org_1 are found in org_1 compared with org_2
            gene, ortho_gene = getGenesPair(line)
            bitscore = getBitScore(line)
            if ortho_gene in gene_pairs:
                if gene == gene_pairs[ortho_gene ]:
                    df.loc[index]=[gene, ortho_gene, bitscore]
                    index += 1
        lines_counter += 1
    source_file2.close()
    df.to_csv(output_path + "/reciprocal_hits_" + org_1 + "_vs_" +org_2 + ".csv", index=False) #pushing data frame to csv file


def getOrganismsName(blast_out_line):
        """
        :param blast_out_line: first line of BLAST output file comparing org1 with org2
        :return (organism_1, organism_2): names of two organisms compared
        """
        res = re.split('\s+', blast_out_line)
        org_1 = re.split('_cds_', res[0][4:-2])[0]
        org_2 = re.split('_cds_', res[1][4:-2])[0]
        return (org_1, org_2)


def getGenesPair(blast_out_line):
    """
    :param blast_out_line: line of BLAST output file comparing org1 with org2
    :return (gene, ortho_gene): pair of orthologous genes.
    """
    res = re.split('\s+', blast_out_line)
    gene = re.split('_cds_',res[0][4:-2])[1]
    ortho_gene = re.split('_cds_',res[1][4:-2])[1]
    return (gene, ortho_gene)

def getBitScore(blast_out_line):
    """
    :param blast_out_line: line of BLAST output file comparing org1 with org2
    :return bit_score of pair of genes compared in the line
    """
    res = re.split('\s+', blast_out_line)
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
        logger = logging.getLogger('main')
        from sys import argv
        if len(argv) < 4:  #
            logger.error('Usage: python ' + argv[0] + ' <arg1> <arg2> <arg3>')
            exit()
        else:
            reciprocalHits(*argv[1:])
