# -*- coding: utf-8 -*-
"""
Running example:
script_name.py /Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/concatenated_all_reciprocal_hits.txt /Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/07_putative_table/putative_orthologs_table.txt /Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/08_mcl_input_files

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
    ** output dir - path to a folder that will contain input files for mcl.

Output:
    Files that are ready to be sent to mcl.

"""

import os

def load_reciprocal_hits_to_dictionary(all_reciprocal_hits_path, group_name_to_pair_combinations, delimiter):
    gene_pair_to_score = {}
    with open(all_reciprocal_hits_path) as f:
        for line in f:
            line_tokens = line.rstrip().split(delimiter)
            if 'bitscore' in line:
                # new reciprocal hits file starts
                continue
            pair = tuple(line_tokens[:2])
            for group_name in group_name_to_pair_combinations:
                pair_combinations = group_name_to_pair_combinations[group_name]
                if pair in pair_combinations or pair[::-1] in pair_combinations:
                    score = line_tokens[2]
                    gene_pair_to_score[pair] = score
                    break
    return gene_pair_to_score


def generate_text_to_mcl_input_file(gene_pair_to_score_dict, gene_pairs):
    first_non_existing_pair = True
    text_to_mcl_file = ''
    logger.debug(f'gene_pair_to_score_dict content is:\n{gene_pair_to_score_dict}')
    for pair in gene_pairs:
        # get gene1,gene2 score and if it's not there try gene2,gene1
        score = gene_pair_to_score_dict.get(pair, gene_pair_to_score_dict.get(pair[::-1]))
        if not score:
            # both pairs are not in the dictionary!
            # E.g., say we found these pairs a1-b1 (A:B), a1-c2 (A:C), b1-c1 (B:C), c1-a2 (C:A)
            # a1-c1 was discarded earlier because it is not *best* hit (best hit was with paralogs, a2 and c2, respectively)
            # Thus, it will not appear in the concatenated best hits.
            if first_non_existing_pair:
                logger.fatal(f'Pair {pair} does not exist in gene_pair_to_score_dict!! Setting score to 1\n'
                             f'Notice that there might be some more non existing pairs (turn on debug mode to get the full log).')
                first_non_existing_pair = False
            else:
                # avoid flooding the logs
                logger.debug(f'Pair {pair} does not exist in gene_pair_to_score_dict!! Setting score to 1')
            score = '1'

        text_to_mcl_file += f'{pair[0]}\t{pair[1]}\t{score}\n'
    return text_to_mcl_file


def prepare_files_for_mcl(all_reciprocal_hits_path, putative_orthologs_path, start, end, output_path, delimiter):
    from itertools import combinations
    group_name_to_combinations = {}
    with open(putative_orthologs_path) as f:
        line_number = 0
        f.readline()  # skip header
        for line in f:
            line_number += 1
            if line_number < start:
                continue
            if line_number > end:
                break
            # parse each table row
            line_tokens = line.rstrip().split(delimiter)
            group_name = line_tokens[0]
            og_members = line_tokens[1:]

            # remove empty tokens
            og_members = [member for member in og_members if member != '']

            # get all pair combinations as a set of tuples
            group_name_to_combinations[group_name] = set(combinations(og_members, 2))

    logger.info('Loading relevant reciprocal hit scores to dictionary...')
    gene_pair_to_score_dict = load_reciprocal_hits_to_dictionary(all_reciprocal_hits_path, group_name_to_combinations, delimiter)
    logger.info(f'All relevant reciprocal hits were loaded successfully. Number of relevant gene pairs is {len(gene_pair_to_score_dict)}.')

    for group_name in group_name_to_combinations:
        # create input file for mcl
        text_to_mcl_file = generate_text_to_mcl_input_file(gene_pair_to_score_dict, group_name_to_combinations[group_name])
        mcl_file_path = os.path.join(output_path, group_name + '.mcl_input')

        with open(mcl_file_path, 'w') as mcl_f:
            mcl_f.write(text_to_mcl_file)


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('all_reciprocal_hits_path', help='path to a file with all the reciprocal hits files concatenated')
    parser.add_argument('putative_orthologs_path', help='path to a file with the putative orthologs sets')
    parser.add_argument('start', help='first ortholog set to prepare', type=int)
    parser.add_argument('end', help='last ortholog set to prepare', type=int)
    parser.add_argument('output_path', help='a folder in which the input files for mcl will be written')
    parser.add_argument('--delimiter', help='delimiter for the input and output files', default=',')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    prepare_files_for_mcl(args.all_reciprocal_hits_path, args.putative_orthologs_path,
                          args.start, args.end, args.output_path, args.delimiter)