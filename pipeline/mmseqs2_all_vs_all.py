def blast_all_vs_all(query_dna_db, query_aa_db, target_dna_db, target_aa_db, aln_db, aln_offsetted_db, tmp_dir, m8_outfile):
    '''
    input:  mmseqs2 DBs
    output: query_vs_reference "mmseqs2 search" results file
    '''
    import subprocess

    cmd = f'mmseqs search {query_aa_db} {target_aa_db} {aln_db} {tmp_dir} --max-seqs 1'
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)

    cmd = f'mmseqs offsetalignment {query_aa_db} {query_dna_db} {target_aa_db} {target_dna_db} {aln_db} {aln_offsetted_db}'
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)

    cmd = f'mmseqs convertalis {query_aa_db} {target_aa_db}  {aln_offsetted_db}  {m8_outfile}'
    logger.info(f'Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('query_dna_db', help='path to mmseqs2 QUERY dna DB')
        parser.add_argument('query_aa_db', help='path to mmseqs2 QUERY aa DB')

        parser.add_argument('target_dna_db', help='path to mmseqs2 TARGET dna DB')
        parser.add_argument('target_aa_db', help='path to mmseqs2 TARGET aa DB')

        parser.add_argument('aln_db', help='path to mmseqs2 alignment DB')
        parser.add_argument('aln_offsetted_db', help='path to mmseqs2 offsetted alignment DB')

        parser.add_argument('tmp_dir', help='path to mmseqs2 intermediate tmp files')
        parser.add_argument('output_path', help='path to which the results will be written (blast m8 format)')

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        blast_all_vs_all(args.query_dna_db, args.query_aa_db, args.target_dna_db, args.target_aa_db,
                         args.aln_db, args.aln_offsetted_db, args.tmp_dir, args.output_path)
