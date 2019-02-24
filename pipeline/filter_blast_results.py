import pandas as pd


def filter_blast(query_vs_reference, output_path, precent_identity_cutoff, e_value_cutoff, delimiter):
    '''
    input:  path to file of blast results
            desiered cutoff values
    output: file with filtered results
    '''
    header = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    df = pd.read_csv(query_vs_reference, header=None, sep=delimiter, names=header.split())
    result = df[(df.pident >= precent_identity_cutoff) & (df.evalue <= e_value_cutoff)]
    result.to_csv(output_path, sep=delimiter, index=False)


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('blast_result', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('--identity_cutoff', help='path to translated sequences', type=float)
    parser.add_argument('--evaule_cutoff', help='path to translated sequences', type=float)
    parser.add_argument('--delimiter', help='orthologs table delimiter', default='\t')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    filter_blast(args.blast_result, args.output_path, args.identity_cutoff, args.evaule_cutoff, args.delimiter)

