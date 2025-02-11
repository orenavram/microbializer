import subprocess
import argparse
from pathlib import Path
import shutil
import sys
from Bio import SeqIO
from ete3 import Tree, TreeStyle

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import add_default_step_args, str_to_bool, run_step
from auxiliaries import consts


def extract_msa_dimensions(msa_path):
    data = list(SeqIO.parse(msa_path, 'fasta'))
    n_seq = len(data)
    n_loci = len(data[0])
    return n_seq, n_loci


def RAxML_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                      bootstrap, seed):
    final_tree_name = phylogenetic_tree_path.name

    final_tree_path = tmp_folder / f'RAxML_result.{final_tree_name}'
    cmd = f'raxmlHPC-PTHREADS-SSE3 -m PROTGAMMAILG -p {seed} -s {msa_path} -n {final_tree_name} -w {tmp_folder} -T {num_of_cpus}'

    if outgroup:
        logger.info(f'The following outgroup was provided: {outgroup}')
        cmd += ' ' + f'-o {outgroup}'
    if bootstrap:
        logger.info(f'{consts.NUMBER_OF_RAXML_BOOTSTRAP_ITERATIONS} bootstrap iterations are going to be done')
        final_tree_path = tmp_folder / f'RAxML_bipartitions.{final_tree_name}'
        cmd += f' -f a -x {seed} -N {consts.NUMBER_OF_RAXML_BOOTSTRAP_ITERATIONS}'

    logger.info(f'Reconstructing species phylogeny with RAxML. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True, check=True)

    if final_tree_path.exists():
        logger.info(f'Copying result {final_tree_path} to {phylogenetic_tree_path}')
        shutil.copy(final_tree_path, phylogenetic_tree_path)
    else:
        logger.fatal(f'TREE WAS NOT GENERATED!!')

    # update info file regarding duplicated sequences reduction
    raxml_info_output_path = tmp_folder / f'RAxML_info.{final_tree_name}'
    if raxml_info_output_path.exists():
        with open(raxml_info_output_path) as f:
            info = f.read()
        info = info.replace(f' {msa_path.parent}/', ': ').replace('reduced',
                                                                  'reduced (at the same folder of the reconstructed core genome)')
        with open(raxml_info_output_path, 'w') as f:
            f.write(info)


def iqtree_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                       bootstrap, seed):
    search_prefix = tmp_folder / 'species_tree'
    tree_output = search_prefix.with_suffix('.treefile')

    cmd = f"iqtree -s {msa_path} -m WAG+G -seed {seed} -pre {search_prefix} -T {num_of_cpus}"
    if outgroup:
        logger.info(f'The following outgroup was provided: {outgroup}')
        cmd += ' ' + f'-o {outgroup}'
    if bootstrap:
        logger.info(f'{consts.NUMBER_OF_IQTREE_BOOTSTRAP_ITERATIONS} bootstrap iterations are going to be done')
        cmd += f' -B {consts.NUMBER_OF_IQTREE_BOOTSTRAP_ITERATIONS}'

    logger.info(f'Reconstructing species phylogeny with IQTree. Executed command is: {cmd}')
    subprocess.run(cmd, shell=True, check=True)
    logger.info(f'IQTree finished successfully. The tree was saved to {tree_output}')

    if tree_output.exists():
        logger.info(f'Copying result {tree_output} to {phylogenetic_tree_path}')
        shutil.copy(tree_output, phylogenetic_tree_path)
    else:
        logger.fatal(f'TREE WAS NOT GENERATED!!')


def fasttree_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                         bootstrap, seed):
    cmd = f"FastTree {msa_path} > {phylogenetic_tree_path}"

    logger.info(f'Reconstructing species phylogeny with FastTree. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True, check=True)


def draw_tree(logger, phylogenetic_tree_path, bootstrap):
    with open(phylogenetic_tree_path, "r") as f:
        tree = Tree(f.readline().strip())

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True

    if bootstrap:
        ts.show_branch_support = True

    tree_image_png_path = phylogenetic_tree_path.with_suffix('.png')
    tree_image_svg_path = phylogenetic_tree_path.with_suffix('.svg')

    logger.info('Drawing the phylogenetic tree. The tree will be saved as PNG and SVG files in the same folder as '
                'the tree file.')
    tree.render(str(tree_image_png_path), tree_style=ts)
    tree.render(str(tree_image_svg_path), tree_style=ts)
    logger.info(f'The tree was saved as PNG and SVG files in the same folder as the tree file.')


def generate_phylogenetic_tree(logger, msa_path, phylogenetic_tree_path, tmp_folder, seed, tree_search_software,
                               outgroup,
                               bootstrap, num_of_cpus):
    # optionally use to decide which program to use
    # n_seq, n_loci = extract_msa_dimensions(msa_path)

    if not phylogenetic_tree_path.exists():  # Might already exist if the job restarted
        if tree_search_software == 'raxml':
            RAxML_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                              bootstrap, seed)
        elif tree_search_software == 'iqtree':
            iqtree_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                               bootstrap, seed)
        elif tree_search_software == 'fasttree':
            fasttree_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                                 bootstrap, seed)

    draw_tree(logger, phylogenetic_tree_path, bootstrap)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('msa_path', type=Path, help='path to a multiple sequence alignment file')
    parser.add_argument('phylogenetic_raw_tree_path', type=Path,
                        help='path to an output file in which the phylogenetic tree will be written')
    parser.add_argument('tmp_path', type=Path, help='path to a tmp folder')
    parser.add_argument('--seed', help='RaxML seed parameter', default=12345)
    parser.add_argument('--tree_search_software', default='iqtree',
                        help='Tree search software to perform phylogenetic tree search. Use iqtree/raxml/fasttree')
    parser.add_argument('--outgroup', help='outgroup for the tree search')
    parser.add_argument('--bootstrap', type=str_to_bool,
                        help='whether or not to apply bootstrap procedure over the reconstructed species tree.')
    parser.add_argument('--cpu', default='1',
                        help='How many CPUs will be used? (for running in parallel mode). For further details see:\nhttps://support.nesi.org.nz/hc/en-gb/articles/115001854444-RAxML#parallel-versions')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, generate_phylogenetic_tree, args.msa_path, args.phylogenetic_raw_tree_path, args.tmp_path, args.seed,
             args.tree_search_software, args.outgroup, args.bootstrap, args.cpu)
