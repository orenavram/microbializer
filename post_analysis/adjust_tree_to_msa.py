import os
import subprocess
import re
from Bio import Phylo
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')



def prune_trees(msa_path, tree_to_prune_path, tmp_dir, output_path):

    with open(tree_to_prune_path) as f:
        match = re.search('\)\d+:',f.read())
        if match:
            logger.error(f'A bootstrap value was detected ->  {match.group()}')
            raise TypeError('Cannot prune tree with bootstrap values. Please provide a tree without them.')

    msa_name = os.path.split(msa_path)[1]
    logger.debug(f'Creating tree for {msa_path}')

    # get list of taxa in the full tree
    tree_taxa = [node.name for node in Phylo.read(tree_to_prune_path, 'newick').get_terminals()]

    # get list of taxa in msa
    with open(msa_path) as f:
        msa = f.read()
    msa_taxa = re.findall(r'>(\S+)\r?\n', msa)

    list_of_taxa_to_remove = [taxon for taxon in tree_taxa if taxon not in msa_taxa]

    # list of names to prune
    list_of_taxa_to_remove_path = f'{tmp_dir}/{msa_name}.txt'
    with open(list_of_taxa_to_remove_path, 'w') as f:
        f.write('\n'.join(list_of_taxa_to_remove))

    cmd_for_prunning = f'/groups/pupko/orenavr2/src/removeTaxa {tree_to_prune_path} {list_of_taxa_to_remove_path} {output_path}'
    subprocess.run(cmd_for_prunning, shell=True)

    logger.info(f'Tree was prunned successfully to {output_path}')



if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('msa_path',
                            help='A path to an MSA file to which the tree should be adjusted',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('tree_to_prune_path',
                            help='A path to a species tree that contains (at least) all the species in the input MSA',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('tmp_dir',
                            help='A path to a folder in which a txt file (with the same name as the msa_file) will be'
                                 'created. All the species that do not appear in the msa (and thus will be removed) '
                                 'will be written to the file that was created',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_path',
                            help='A path to a file in which the prunned tree will be written',
                            type=lambda path: path if os.path.exists(os.path.split(path)[0]) else parser.error(
                                f'output folder {os.path.split(path)[0]} does not exist!'))
        args = parser.parse_args()

        prune_trees(args.msa_path, args.tree_to_prune_path, args.tmp_dir, args.output_path)
