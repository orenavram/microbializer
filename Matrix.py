import pandas as pd
import csv


def generateOrthologousMatrix(matrix_path, reciprocal_hits_path):
    """
    :param matrix_path: path of OMA (Orthologous Matrix) file to be revised (if not exist- pass with desired file name)
    :param  reciprocal_hits_path: path of reciprocal_hits file
    updating global OMA csv file basing on genes & Orthologous genes of organisms in reciprocal hits file
    """
    OMA = editCsvFile(matrix_path)  #generating 'editable' data frame from matrix file
    reciprocal_hits = readCsvFile(reciprocal_hits_path)    #generating 'read only' data frame from Symetric blast file
    org_1, org_2 =getOrganismsName(reciprocal_hits)   #getting names of organisms compared in Symetric blast file
    gene_index = -1
    ortho_gene_index = -1
    for index, row in reciprocal_hits.iterrows():
        gene,ortho_gene = getGenePair(row, org_1, org_2)
        addGenesToMat(OMA, gene, ortho_gene, org_1, org_2)
    OMA.to_csv(matrix_path, index=False) #pushing revised OMA data frame to global csv file


def editCsvFile(matrix_path):
    """
    :param matrix_path: path of OMA (Orthologous Matrix) file to be revised
    :return: data frame of OMA  to be edited.
    """
    try: #matrix file (csv) already exist.
        df = pd.read_csv(matrix_path, na_filter=False)
    except Exception: #matrix file doesn't exist - creating new csv file.
        with open (matrix_path, 'w', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(["GENE"])   #adding the coulmn "GENE" to OMA
    return pd.read_csv(matrix_path, na_filter=False)

def readCsvFile(reciprocal_hits_path):
    """
    :param reciprocal_hits_path: path of reciprocal hits file to be read
    :return: data frame of reciprocal hits to be read and analyzed.
    """
    return pd.read_csv(reciprocal_hits_path)

def getOrganismsName(reciprocal_hits):
    """
    :param reciprocal_hits: data frame of reciprocal hits.
    :return (organism_1, organism_2): names of two organisms compared in the data frame
    """
    orgs_list = list(reciprocal_hits.columns.values);
    return (orgs_list[0], orgs_list[1])

def getGenePair(row, organism_1, organism_2):
    """
    :param row: row in reciprocal_hits data frame
    :param organism_1, organism_2: names of two organisms compared in the data frame
    :return (gene, ortholog_gene): pair of orthologous genes.
    """
    return (row[organism_1], row[organism_2])

def addGenesToMat(OMA_df, gene, ortho_gene, org_1, org_2):
    """
    :param OMA_df: data frame of OMA to be revised
    :param gene, ortho_gene: pair of orthologous genes.
    :param org_1, org_2: names of two organisms containig gene and ortho_gene accordingly
    adding genes & Orthologous genes of organisms to OMA data frame

    """
    gene_index= -1
    ortho_gene_index = -1
    if gene in OMA_df.GENE.values: #checking if gene is found on "GENE" coulmn & updating orgs coulmns accordingly
        gene_index = OMA_df[OMA_df.GENE == gene].index[0]
        OMA_df.at[gene_index, org_2] = ortho_gene
        OMA_df.at[gene_index, org_1] = gene
    elif ortho_gene in  OMA_df.GENE.values: #checking if ortho_gene is found on "GENE" coulmnu & updating orgs coulmns accordingly
        ortho_gene_index = OMA_df[OMA_df.GENE == ortho_gene].index[0]
        OMA_df.at[ortho_gene_index, org_1] = gene
        OMA_df.at[ortho_gene_index, org_2] = ortho_gene
    else: #adding gene to "GENE" coulmn & updating orgs coulmns accordingly
        OMA_df.at[-1, "GENE"] = gene
        OMA_df.at[-1, org_2] = ortho_gene
        OMA_df.at[-1, org_1] = gene
        OMA_df.index += 1



if __name__ == '__main__':
        import logging
        logger = logging.getLogger('main')
        from sys import argv
        if len(argv) < 3:  #
            logger.error('Usage: python ' + argv[0] + ' <arg1> <arg2>')
            exit()
        else:
            generateOrthologousMatrix(*argv[1:])