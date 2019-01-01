def verify(input_file, output_file, clustering_criterion):
    with open(input_file) as f:
        i = 0
        for line in f:
            # each line is a single cluster
            i += 1
            if i > clustering_criterion:
                logger.info(f'{input_file} has at least {i} clusters, thus not relevant.')
                return False
        if i == 0:
            raise ValueError(f'{input_file} is empty! There\'s a bug in the previous step!')
    import os
    os.rename(input_file, output_file)
    return True



if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_file', help='path to an MCL analysis file')
        parser.add_argument('output_file', help='path to which the MCL analysis will be moved if clustering criterion was met')
        parser.add_argument('--clustering-criterion', help='maximal number of clusters allowed', type=int, default=1)
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        verify(args.input_file, args.output_file, args.clustering_criterion)








