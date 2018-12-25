'''
script_name.py /Users/Oren/Dropbox/Projects/microbializer/mock_output/04_blast_filtered/GCF_000009505.1_ASM950v1_cds_from_genomic_vs_GCF_000195995.1_ASM19599v1_cds_from_genomic.04_blast_filtered /Users/Oren/Dropbox/Projects/microbializer/mock_output/04_blast_filtered/GCF_000195995.1_ASM19599v1_cds_from_genomic_vs_GCF_000009505.1_ASM950v1_cds_from_genomic.04_blast_filtered /Users/Oren/Dropbox/Projects/microbializer/mock_output/05_reciprocal_hits/GCF_000009505.1_ASM950v1_cds_from_genomic_vs_GCF_000195995.1_ASM19599v1_cds_from_genomic.05_reciprocal_hits
'''

import os


def parse_blast_results_to_dictionary(blast_out):
    query_to_hit_and_bitscore_dict = {}
    with open(blast_out) as f:
        line = f.readline()  # skip header
        assert 'bitscore' in line, f'A non header row was skipped!\n{line}'
        for line in f:
            line_tokens = line.rstrip().split()
            query, hit, bitscore = line_tokens[0], line_tokens[1], line_tokens[-1]
            query_to_hit_and_bitscore_dict[query] = (hit, bitscore)
    return query_to_hit_and_bitscore_dict


def find_reciprocal_hits(blast_out1, blast_out2, output_path, delimiter):
    output_file_name = os.path.split(output_path)[1]
    strain1, strain2 = os.path.splitext(output_file_name)[0].split('_vs_')
    result = delimiter.join([strain1, strain2, 'bitscore']) + '\n'

    query_to_hit_and_bitscore_dict_1 = parse_blast_results_to_dictionary(blast_out1)
    query_to_hit_and_bitscore_dict_2 = parse_blast_results_to_dictionary(blast_out2)

    for query1 in query_to_hit_and_bitscore_dict_1:
        hit1, bitscore1 = query_to_hit_and_bitscore_dict_1[query1]
        # swap query and hit for the other blast file
        query2, hit2 = hit1, query1
        query2_val = query_to_hit_and_bitscore_dict_2.get(query2, (None, None))
        if query2_val[0] == hit2:
            # a reciprocal hit was found!
            result += delimiter.join([query1, query2, bitscore1]) + '\n'

            # sanity check that reciprocal bitscores are identical
            bitscore2 = query2_val[1]
            if eval(bitscore1) != eval(bitscore2):
                logger.warning(f'bitscores of {query1} and {query2} are different!!\n{bitscore1} in {blast_out1}\n{bitscore2} in {blast_out2}')

            #TODO: remove if causes crash!
            assert eval(bitscore1) == eval(bitscore2), f'bitscores of {query1} and {query2} are different!!\n{bitscore1} in {blast_out1}\n{bitscore2} in {blast_out2}'

    with open(output_path, 'w') as f:
        f.write(result)
    logger.debug(f'{output_path} was written seccessfully')


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('blast_result1', help='path to blast result file of seq1 vs seq2')
    parser.add_argument('blast_result2', help='path to blast result file of seq2 vs seq1')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=',')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    find_reciprocal_hits(args.blast_result1, args.blast_result2, args.output_path, args.delimiter)
