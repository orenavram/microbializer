def search_all_vs_all(aa_db1, aa_db2, aln_offsetted_db, tmp_dir, m8_outfile):
    '''
    input:  mmseqs2 DBs
    output: query_vs_reference "mmseqs2 search" results file
    '''
    import os
    import subprocess

    i = 1
    while not os.path.exists(m8_outfile):
        # when the data set is very big some files are not generated because of the heavy load
        # so we need to make sure they will be generated!
        logger.info(f'Iteration #{i}: Executing mmseqs2 search for {aa_db1} and {aa_db2}')
        i += 1

        # control verbosity level by -v [3] param ; verbosity levels: 0=nothing, 1: +errors, 2: +warnings, 3: +info
        v = 2

        cmd = f'mmseqs rbh {aa_db1} {aa_db2} {aln_offsetted_db} {tmp_dir} -v {v}'

        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)

        cmd = f'mmseqs convertalis {aa_db1} {aa_db2}  {aln_offsetted_db}  {m8_outfile} -v {v}'
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)

    intermediate_files_prefix = os.path.splitext(aln_offsetted_db)[0]
    intermediate_files = [f'{intermediate_files_prefix}{suffix}' for suffix in ['.alnOffsettedDB', '.alnOffsettedDB.dbtype', '.alnOffsettedDB.index']]
    # each pair generates 3 intermidiate files! lot's of junk once finished
    for file in intermediate_files:
        os.remove(file)


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('aa_db1', help='path to an aa DB')
        parser.add_argument('aa_db2', help='path to another aa DB')
        parser.add_argument('aln_offsetted_db', help='path to mmseqs2 offsetted alignment DB')
        parser.add_argument('tmp_dir', help='a path to write mmseqs internal files')
        parser.add_argument('output_path', help='path to which the results will be written (blast m8 format)')

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        search_all_vs_all(args.aa_db1, args.aa_db2, args.aln_offsetted_db,
                          args.tmp_dir, args.output_path)
