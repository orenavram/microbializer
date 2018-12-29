def get_verified_clusters_set(verified_clusters_path):
    import os
    return set([os.path.splitext(file)[0] for file in os.listdir(verified_clusters_path)])


def finalize_table(putative_orthologs_path, verified_clusters_path, finalized_table_path, table_header_path, delimiter):
    verified_clusters_set = get_verified_clusters_set(verified_clusters_path)
    logger.info(f'verified_clusters_set:\n{verified_clusters_set}')
    with open(putative_orthologs_path) as f:
        # OG_name,GCF_000008105,GCF_000006945,GCF_000195995,GCF_000007545,GCF_000009505
        putative_orthologs_table_header = f.readline().rstrip()
        # remove "OG_name," from header
        index_of_first_delimiter = putative_orthologs_table_header.index(delimiter)
        final_orthologs_table_header = putative_orthologs_table_header[index_of_first_delimiter + 1:]
        finalized_table_str = final_orthologs_table_header + '\n'
        for line in f:
            first_delimiter_index = line.index(delimiter)
            OG_name = line[:first_delimiter_index]
            group = line.rstrip()[first_delimiter_index+1:]
            if OG_name in verified_clusters_set:
                logger.debug(f'Adding {OG_name} to the final table')
                verified_clusters_set.discard(OG_name)
                finalized_table_str += group + '\n'

    with open(table_header_path, 'w') as f:
        f.write(final_orthologs_table_header)

    with open(finalized_table_path, 'w') as f:
        f.write(finalized_table_str)



if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('putative_orthologs_path', help='path to a file with the putative orthologs sets')
        parser.add_argument('verified_clusters_path', help='path to a directory with the verified clusters')
        parser.add_argument('finalized_table_path',
                            help='path to an output file in which the final table will be written')
        parser.add_argument('table_header_path',
                            help='path to an output file in which ONLY THE HEADER of the final table will be written')
        parser.add_argument('--delimiter', help='delimiter for the putative orthologs table', default=',')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        finalize_table(args.putative_orthologs_path, args.verified_clusters_path, args.finalized_table_path,
                       args.table_header_path, args.delimiter)








