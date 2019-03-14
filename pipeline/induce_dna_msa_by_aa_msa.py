"""

"""

def induce_sequence(aa_seq, dna_seq):
    result = ''
    for i in range(len(aa_seq)):
        if aa_seq[i] == '-':
            result += '-'*3
        else:
            result += dna_seq[i*3:(i+1)*3]

    return result


def induce_msa(aa_msa_path, dna_ms_path, output_path):
    from extract_core_genome import load_orthologs_group_to_dict
    og_name_to_aa = load_orthologs_group_to_dict(aa_msa_path)
    og_name_to_dna = load_orthologs_group_to_dict(dna_ms_path)

    result = ''
    with open(dna_ms_path) as f:
        for line in f:
            if line.startswith('>'):
                og_name = line.lstrip('>').rstrip()
                induced_dna_sequence = induce_sequence(og_name_to_aa[og_name], og_name_to_dna[og_name])
                result += f'>{og_name}\n{induced_dna_sequence}\n'

    with open(output_path, 'w') as f:
        f.write(result)


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('aa_msa_path', help='path to a file with aligned proteins')
    parser.add_argument('dna_ms_path', help='path to a file with unaligned dna sequences')
    parser.add_argument('output_path', help='path to a file in which the induced dna alignment will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    induce_msa(args.aa_msa_path, args.dna_ms_path, args.output_path)


