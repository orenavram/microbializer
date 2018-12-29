import matplotlib
matplotlib.use('Agg')
from Bio import Phylo
import pylab as pb


def generatePhylogeneticTree(msa_path, phylogenetic_tree_path, seed, model):
    import os
    wd, phylogenetic_tree_name = os.path.split(phylogenetic_tree_path)

    cmd = f'raxmlHPC -s {msa_path} -w {wd} -n {phylogenetic_tree_name} -p {seed} -m {model}'
    import subprocess
    subprocess.run(cmd, shell=True)

    raxml_best_tree_output_path = os.path.join(wd, f'RAxML_bestTree.{phylogenetic_tree_name}')
    if os.path.exists(raxml_best_tree_output_path):
        os.rename(raxml_best_tree_output_path, phylogenetic_tree_path)
    else: # Tree was NOT generated!
        #TODO: make sure that the script of ploting the tree is aware to that there is no tree.
        pass


    raxml_info_output_path = os.path.join(wd, f'RAxML_info.{phylogenetic_tree_name}')
    if os.path.exists(raxml_info_output_path):
        with open(raxml_info_output_path) as f:
            info = f.read()
        info = info.replace(f' {os.path.split(msa_path)[0]}/', ': ').replace('reduced', 'reduced (at the same folder of the reconstructed core genome)')
        with open(raxml_info_output_path, 'w') as f:
            f.write(info)


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('msa_path', help='path to a multiple sequence alignment file')
        parser.add_argument('phylogenetic_tree_path', help='path to an output file in which the phylogenetic tree will be written')
        parser.add_argument('--seed',
                            help='RaxML seed parameter',
                            default=12345)
        parser.add_argument('--model',
                            help='RaxML reconstruction model parameter',
                            default='GTRGAMMAI')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        generatePhylogeneticTree(args.msa_path, args.phylogenetic_tree_path,
                            args.seed, args.model)




def generateTreeFigure(tree_file_directory, output_file_path='tree_figure'):
    """Generating figure based on phylogenetic tree file.

    Parameters
    ----------
    tree_file_directory : string
       Path of file containing newick format tree.
    output_file_path: string
        Path of figure output file.

    Returns
    -------
    None
        generating tree figure.
    """
    tree = Phylo.read(tree_file_directory, 'newick')
    Phylo.draw(tree, do_show=False)
    pb.savefig(output_file_path)
