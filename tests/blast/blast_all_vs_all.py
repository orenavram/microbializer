from sys import argv
import argparse
import logging
import subprocess


def blast_all_vs_all(logger, program, query_seq, reference_DB, output_path):
    """
    input:  type of blast to run(blastn/blastp)
            fasta file of query genome
            blast DB file
    output: query_vs_reference blast results file
    """
    cmd = f'{program} -query {query_seq} -db {reference_DB} -out {output_path} -max_target_seqs 1 -max_hsps 1 -outfmt 6'
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    parser = argparse.ArgumentParser()
    parser.add_argument('query_fasta', help='path to query fasta file')
    parser.add_argument('subject_db', help='path to subject blast DB file that will be created')
    parser.add_argument('output_path', help='path to output file - will contain the blast results')
    parser.add_argument('--blast_type', help='blast type to run (blastn/blastp)',
                        default='blastn')  # TODO: why not tblastx?
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    blast_all_vs_all(logger, args.blast_type, args.query_fasta, args.subject_db, args.output_path)
