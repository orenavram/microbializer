import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def fix_column(column, strain2sequence, legal_chars='ACGT-'):
    col = [strain2sequence[strain][column] for strain in strain2sequence
           if strain2sequence[strain][column] in legal_chars]

    major_allele = max(set(col), key=col.count)
    for strain in strain2sequence:
        if strain2sequence[strain][column] not in legal_chars:
            logger.info(f'Replacing column #{column} of {strain} from {strain2sequence[strain][column]} to {major_allele}')
            strain2sequence[strain] = strain2sequence[strain][:column] + major_allele + strain2sequence[strain][column+1:]


def main(msa_path, tree_path, control_path, output_path,
         MPReconstruct_script_path='/groups/pupko/orenavr2/pupkoSVN/programs/MPreconstruct/MPreconstruct'):

    # make sure that an appropriate gcc module is loaded. E.g., gcc/gcc-6.2.0
    import subprocess

    # prepare a control file for the c++ code:
    with open(control_path, 'w') as f:
        f.write(f'''### parameters for MPreconstruct

## in and out files
_treefile {tree_path}
_seqfile {msa_path}
_logfile {os.path.splitext(control_path)[0]}.log
_outfile {os.path.splitext(control_path)[0]}.out  

## types are: nuc, amino, threeState, integer
_alphabetType nuc

## alphabet: a, c, g, t, and -
_alphabetSize 5

## types are: file,fitch,diff,diffSquare
_costMatrixType fitch''')

    # what to execute in the shell:
    cmd = f'{MPReconstruct_script_path} {control_path} {output_path}'

    logger.info(f'Computing homoplasy. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('msa_path',
                            help='A path to an MSA file to compute homoplasy',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('tree_path',
                            help='A path to a(n adjusted) phylogenetic tree',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('control_path',
                            help='A path to control file in which the MPReconstruct parameters will be written to',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        parser.add_argument('output_path',
                            help='A path in which the computed homoplasy will be written to',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        main(args.msa_path, args.tree_path, args.control_path, args.output_path)


