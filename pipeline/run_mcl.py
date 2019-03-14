def mcl(input_file, output_file):
    # --abc for a columns format, i.e., item1\item2\tscore
    import subprocess
    cmd = f'mcl {input_file} --abc -o {output_file}'
    logger.info(f'Starting MCL. Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_file', help='path to an MCL input file')
        parser.add_argument('output_file', help='path to which the MCL analysis will be written')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        mcl(args.input_file, args.output_file)





