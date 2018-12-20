import subprocess

def create_blast_DB(reference_seq, dbtype, output_path):
    '''
    input:  sequnce to base the DB on
            DB type (nucl/prot)
            path to output file
    output: balst DB based on the reference sequence
    '''
    cmd = f'makeblastdb -in {reference_seq} -out {output_path} -dbtype {dbtype}'
    logger.info(f'Starting makeblastdb. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    import logging
    logger = logging.getLogger('main')
    from sys import argv
    logger.info(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('subject_fasta', help='path to subject fasta file')
    parser.add_argument('subject_db', help='path to subject blast DB file that will be created')
    parser.add_argument('--dbtype', help='database type for creating blast DB', default='nucl')

    args = parser.parse_args()

    create_blast_DB(args.subject_fasta, args.dbtype, args.subject_db)