def create_mmseq2_DB(fasta_path, output_path, tmp_path, translate, convert2fasta):
    '''
    input:  sequnce to base the DB on
            DB type (nucl/prot)
            path to output file
    output: mmseqs2 DB based on the reference sequence
    '''
    #for more details see: mmseqs createdb -h
    import subprocess

    cmd = f'mmseqs createdb {fasta_path} {tmp_path} --dont-split-seq-by-len'
    logger.info(f'Starting mmseqs createdb. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)

    if translate or convert2fasta:
        cmd = f'mmseqs translatenucs {tmp_path} {tmp_path}_aa'  # translate dna db to aa db
        logger.info(f'Translating dna DB to amino acids. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)

    if convert2fasta:
        cmd = f'mmseqs convert2fasta {tmp_path}_aa {output_path}_aa.fas'  # convert the aa db to a regular fasta
        logger.info(f'Converting amino acids DB to a fasta format. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)



if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fasta', help='path to input fasta file')
    parser.add_argument('output_prefix', help='path prefix for the DB file(s)')
    parser.add_argument('tmp_prefix', help='path prefix for the tmp file(s)')
    #parser.add_argument('--dbtype', help='database type for creating blast DB', default='nucl')
    parser.add_argument('-t', '--translate', help='whether to translate the dna to aa', action='store_true')
    parser.add_argument('-c', '--convert2fasta', help='whether to convert the dbs to fasta', action='store_true')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    create_mmseq2_DB(args.input_fasta, args.output_prefix, args.tmp_prefix, args.translate, args.convert2fasta)