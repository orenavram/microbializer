"""

"""

def reconstruct_msa(sequences_file_path, output_file_path, maxiterate):
    import subprocess
    cmd = f'mafft-linsi --maxiterate {maxiterate} --localpair {sequences_file_path} > {output_file_path}'
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('sequences_file_path', help='path to a file with unaligned sequences')
    parser.add_argument('output_file_path', help='path to a file in which the aligned sequences will be written')
    parser.add_argument('--maxiterate', help='number for MAFFT maxiterate parameter', default=1000, type=int)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    reconstruct_msa(args.sequences_file_path, args.output_file_path, args.maxiterate)


