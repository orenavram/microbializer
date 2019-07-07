import subprocess


# module load prodigal/prodigal-2.6.3
def find_genes(genome, output_path, log_path):
    '''
    input:path to fasta file with prokaryotic genome to be analyzed
    output: protein-coding gene prediction for input genome
    '''
    # segments = list(SeqIO.parse(genome, 'fasta'))
    # length = sum(len(segment) for segment in segments)
    cmd = f'prodigal -i "{genome}"  -d {output_path}'
    logger.info(f'Starting prodigal. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('genome_path', help='path to fasta genome file')
    parser.add_argument('output_path', help='path to output file')
    parser.add_argument('log_path', help='path to translated sequences')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    find_genes(args.genome_path, args.output_path, args.log_path)

