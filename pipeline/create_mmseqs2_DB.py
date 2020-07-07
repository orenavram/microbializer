import os


def too_many_trials(cmd, error_file_path):
    msg = f'Failed to fetch <i>{cmd}</i> command. Could be due to heavy load on our web servers. ' \
          'Please contact us for further assistance.'
    if os.path.exists('/bioseq'):  # remote run
        from sys import path
        path.append('/bioseq/microbializer/auxiliaries')
        from pipeline_auxiliaries import fail
        # get error_log path
        # e.g., from this aa_db1: /bioseq/data/results/microbializer/159375410340094617808216800611/outputs/02_dbs/SAL_BA5690AA_AS.scaffold_aa
        # into this: /bioseq/data/results/microbializer/159375410340094617808216800611/error.txt
        logger.info(f'Writing to error file in {error_file_path}')
        fail(msg, error_file_path)
    raise OSError(msg)


def create_mmseq2_DB(fasta_path, output_path, tmp_path, translate, convert2fasta, verbosity_level):
    '''
    input:  sequnce to base the DB on
            DB type (nucl/prot)
            path to output file
    output: mmseqs2 DB based on the reference sequence
    '''
    #for more details see: mmseqs createdb -h
    import subprocess
    import os
    import time

    i = 1
    while not os.path.exists(tmp_path):
        logger.info(f'Iteration #{i}: createdb. Result should be at {tmp_path}')
        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        cmd = f'mmseqs createdb {fasta_path} {tmp_path} -v {verbosity_level}'# --dont-split-seq-by-len'
        logger.info(f'Starting mmseqs createdb. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)
        i += 1
        if i == 1000:
            error_file_path = f'{os.path.split(output_path)[0]}/../../error.txt'
            too_many_trials('mmseqs createdb', error_file_path)
        time.sleep(1)

    if translate or convert2fasta:
        i = 1
        while not os.path.exists(f'{tmp_path}_aa'):
            logger.info(f'Iteration #{i}: translatenucs. Result should be at {tmp_path}_aa')
            cmd = f'mmseqs translatenucs {tmp_path} {tmp_path}_aa -v {verbosity_level}'  # translate dna db to aa db
            logger.info(f'Translating dna DB to amino acids. Executed command is:\n{cmd}')
            subprocess.run(cmd, shell=True)
            i += 1
            if i == 1000:
                error_file_path = f'{os.path.split(output_path)[0]}/../../error.txt'
                too_many_trials('mmseqs translatenucs', error_file_path)
            time.sleep(1)

    if convert2fasta:  # for step 13: dna to aa
        i = 1
        while not os.path.exists(f'{output_path}_aa'):
            logger.info(f'Iteration #{i}: convert2fasta. Result should be at {output_path}_aa')
            cmd = f'mmseqs convert2fasta {tmp_path}_aa {output_path}_aa.fas -v {verbosity_level}'  # convert the aa db to a regular fasta
            logger.info(f'Converting amino acids DB to a FASTA format. Executed command is:\n{cmd}')
            subprocess.run(cmd, shell=True)
            i += 1
            if i == 1000:
                error_file_path = f'{os.path.split(output_path)[0]}/../../error.txt'
                too_many_trials('mmseqs convert2fasta', error_file_path)
            time.sleep(1)

        intermediate_files = [f'{tmp_path}{suffix}' for suffix in ['', '_aa', '_aa.dbtype', '_aa_h', '_aa_h.dbtype', '_aa_h.index', '_aa.index', '.dbtype', '_h', '_h.dbtype', '_h.index', '.index', '.lookup']]
        # each pair generates 13 intermidiate files!!! lot's of junk once finished
        for file in intermediate_files:
            os.remove(file)



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

    create_mmseq2_DB(args.input_fasta, args.output_prefix, args.tmp_prefix,
                     args.translate, args.convert2fasta, 3 if args.verbose else 1)