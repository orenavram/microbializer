import pandas as pd


def filter_blast(query_vs_reference, output_path, precent_identity_cutoff=80, e_value_cutoff=0.01):
    '''
    input:  path to file of blast results
            desired cutoff values
    output: file with filtered results
    '''
    header = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    df = pd.read_csv(query_vs_reference, header=None, sep='\t', names=header.split())
    result = df[(df.pident >= precent_identity_cutoff) & (df.evalue <= e_value_cutoff)]
    result.to_csv(output_path, sep='\t', index=False)
