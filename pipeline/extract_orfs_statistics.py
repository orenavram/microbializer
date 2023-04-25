def extract_orfs_statistics(orf_path, orfs_count_output_path, orfs_gc_output_path):
    num_of_nucleotides = 0
    num_of_GC = 0
    orfs_count = -1
    sequence = ''
    with open(orf_path) as f:
        for line in f:
            if line.startswith('>'):
                orfs_count += 1
                num_of_nucleotides += len(sequence)
                num_of_GC += sequence.count('G') + sequence.count('C')
                sequence = ''
                if orfs_count % 1000 == 0:
                    logger.debug(f'ORFs count is: {orfs_count}')
            else:
                sequence += line.rstrip().upper()

        # don't forget last record!!
        orfs_count += 1
        num_of_nucleotides += len(sequence)
        num_of_GC += sequence.count('G') + sequence.count('C')

        with open(orfs_count_output_path, 'w') as f:
            f.write(f'{orfs_count}\n')  # [f'{count},{orfs_counts[count]}' for count in orfs_counts]))

        with open(orfs_gc_output_path, 'w') as f:
            f.write(f'{num_of_GC / num_of_nucleotides}\n')


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('orf_path', help='path to fasta file with orfs')
    parser.add_argument('orfs_count_output_path', help='where to write the number of orfs')
    parser.add_argument('orfs_gc_output_path', help='where to write the gc content')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    extract_orfs_statistics(args.orf_path, args.orfs_count_output_path, args.orfs_gc_output_path)
