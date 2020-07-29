import os
import logging

import sys
if os.path.exists('/bioseq'):  # remote run
    sys.path.append('/bioseq/microbializer/auxiliaries')
    from pipeline_auxiliaries import load_header2sequences_dict
else:
    from auxiliaries.pipeline_auxiliaries import load_header2sequences_dict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def fix_column(column, strain2sequence, legal_chars='ACGT-'):
    col = [strain2sequence[strain][column] for strain in strain2sequence
           if strain2sequence[strain][column] in legal_chars]

    major_allele = max(set(col), key=col.count)
    for strain in strain2sequence:
        if strain2sequence[strain][column] not in legal_chars:
            logger.info(f'Replacing column #{column} of {strain} from {strain2sequence[strain][column]} to {major_allele}')
            strain2sequence[strain] = strain2sequence[strain][:column] + major_allele + strain2sequence[strain][column+1:]


def fix_msa(msa_path, output_path):
    strain2sequence, msa_length = load_header2sequences_dict(msa_path, get_length=True, upper_sequence=True)

    for strain in strain2sequence:
        if len(strain2sequence[strain]) != msa_length:
            raise ValueError(f'Illegal MSA. Not all sequences are of the same length. E.g., {strain} sequence length is {len(strain2sequence[strain])} where others are of length {msa_length}')

    for col in range(msa_length):
        fix_column(col, strain2sequence)

    fixed_alignment = ''
    for strain in strain2sequence:
        fixed_alignment += '>' + strain + '\n' + strain2sequence[strain] + '\n'

    with open(output_path, 'w') as f:
        f.write(fixed_alignment)

    logger.info(f'MSA was fixed and stored successfully at {output_path}')


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_msa_path',
                            help='A path to an MSA file that should be fixed',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_msa_path',
                            help='A path to the fixed MSA',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        fix_msa(args.input_msa_path, args.output_msa_path)


