import os


def fna_to_faa(nucleotide_path, protein_path):
    import Bio.Seq
    with open(nucleotide_path) as f:
        result = ''
        previous_protein = ''
        previous_header = f.readline()
        for line in f:
            if line.startswith('>'):
                result += f'{previous_header}{Bio.Seq.translate(previous_protein)}\n'
                previous_protein = ''
                previous_header = line
            else:
                previous_protein += line.rstrip()
        # don't forget to add last record!
        result += f'{previous_header}{Bio.Seq.translate(previous_protein)}\n'

    with open(protein_path, 'w') as f:
        f.write(result)

    logger.info(f'Translated fatsa file was write successfully to: {protein_path}')


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('nucleotide_path',
                        help='A path to a nucleotide fasta file',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('protein_path', help='A path to which the translated dna will be written')
    # type=lambda path: path if os.path.exists(os.path.split(path)[0]) else
    # parser.error(f'output folder {os.path.split(path)[0]} does not exist!'))
    args = parser.parse_args()

    import logging

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    # wait_for_output_folder(os.path.split(args.protein_path)[0])

    fna_to_faa(args.nucleotide_path, args.protein_path)
