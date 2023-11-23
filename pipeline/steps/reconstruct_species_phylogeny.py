import subprocess
from sys import argv
import argparse
import logging
import os
import shutil
import sys
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def extract_msa_dimensions(msa_path):
    data = list(SeqIO.parse(msa_path, 'fasta'))
    n_seq = len(data)
    n_loci = len(data[0])
    return n_seq, n_loci


def RAxML_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                      num_of_bootstrap_iterations, seed):
    final_tree_name = os.path.split(phylogenetic_tree_path)[1]

    final_tree_path = os.path.join(tmp_folder, f'RAxML_result.{final_tree_name}')
    cmd = f'raxmlHPC-PTHREADS-SSE3 -m PROTGAMMAILG -p {seed} -s {msa_path} -n {final_tree_name} -w {tmp_folder} -T {num_of_cpus}'

    if outgroup:
        logger.info(f'The following outgroup was provided: {outgroup}')
        cmd += ' ' + f'-o {outgroup}'
    if num_of_bootstrap_iterations > 0:
        logger.info(f'{num_of_bootstrap_iterations} bootstrap iterations are going to be done')
        final_tree_path = os.path.join(tmp_folder, f'RAxML_bipartitions.{final_tree_name}')
        cmd += f' -f a -x {seed} -N {num_of_bootstrap_iterations}'

    logger.info(f'Reconstructing species phylogeny with RAxML. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)

    if os.path.exists(final_tree_path):
        logger.info(f'Copying result {final_tree_path} to {phylogenetic_tree_path}')
        shutil.copy(final_tree_path, phylogenetic_tree_path)
    else:
        logger.fatal(f'TREE WAS NOT GENERATED!!')

    # update info file regarding duplicated sequences reduction
    raxml_info_output_path = os.path.join(tmp_folder, f'RAxML_info.{final_tree_name}')
    if os.path.exists(raxml_info_output_path):
        with open(raxml_info_output_path) as f:
            info = f.read()
        info = info.replace(f' {os.path.split(msa_path)[0]}/', ': ').replace('reduced',
                                                                             'reduced (at the same folder of the reconstructed core genome)')
        with open(raxml_info_output_path, 'w') as f:
            f.write(info)


def iqtree_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                       num_of_bootstrap_iterations, seed):
    final_tree_name = os.path.split(phylogenetic_tree_path)[1]
    search_prefix = os.path.join(tmp_folder, final_tree_name)

    cmd = f"iqtree -s {msa_path} -m WAG+G -seed {seed} -pre {search_prefix} -T {num_of_cpus}"
    if outgroup:
        logger.info(f'The following outgroup was provided: {outgroup}')
        cmd += ' ' + f'-o {outgroup}'
    if num_of_bootstrap_iterations > 0:
        if num_of_bootstrap_iterations < 1000:
            num_of_bootstrap_iterations = 1000
            logger.info(f'Updating num_of_bootstrap_iterations to 1000 since it is the minimum.')

        logger.info(f'{num_of_bootstrap_iterations} bootstrap iterations are going to be done')
        cmd += f' -B {num_of_bootstrap_iterations}'

    logger.info(f'Reconstructing species phylogeny with IQTree. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)

    final_tree_path = search_prefix + ".treefile"
    if os.path.exists(final_tree_path):
        logger.info(f'Copying result {final_tree_path} to {phylogenetic_tree_path}')
        shutil.copy(final_tree_path, phylogenetic_tree_path)
    else:
        logger.fatal(f'TREE WAS NOT GENERATED!!')


def fasttree_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                         num_of_bootstrap_iterations, seed):
    cmd = f"FastTree {msa_path} > {phylogenetic_tree_path}"

    logger.info(f'Reconstructing species phylogeny with FastTree. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)


def generate_phylogenetic_tree(logger, msa_path, phylogenetic_tree_path, tmp_folder, seed, tree_search_software, outgroup,
                               num_of_bootstrap_iterations, num_of_cpus):
    # optionally use to decide which program to use
    # n_seq, n_loci = extract_msa_dimensions(msa_path)

    if tree_search_software == 'raxml':
        RAxML_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus,outgroup,
                          num_of_bootstrap_iterations, seed)
    elif tree_search_software == 'iqtree':
        iqtree_tree_search(tmp_folder, msa_path, phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                           num_of_bootstrap_iterations, seed)
    elif tree_search_software == 'fasttree':
        fasttree_tree_search(tmp_folder, msa_path,phylogenetic_tree_path, logger, num_of_cpus, outgroup,
                             num_of_bootstrap_iterations, seed)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('msa_path', help='path to a multiple sequence alignment file')
    parser.add_argument('phylogenetic_raw_tree_path',
                        help='path to an output file in which the phylogenetic tree will be written')
    parser.add_argument('tmp_path', help='path to a tmp folder')
    parser.add_argument('--seed',
                        help='RaxML seed parameter',
                        default=12345)
    parser.add_argument('--tree_search_software', default='raxml',
                        help='Tree search software to perform phylogenetic tree search. Use iqtree/raxml/fasttree')
    parser.add_argument('--outgroup', default=None)
    parser.add_argument('--num_of_bootstrap_iterations', type=int, default=0)
    parser.add_argument('--cpu', choices=[str(i) for i in range(1, 29)], default='1',
                        help='How many CPUs will be used? (for running in parallel mode). For further details see:\nhttps://support.nesi.org.nz/hc/en-gb/articles/115001854444-RAxML#parallel-versions')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        generate_phylogenetic_tree(logger, args.msa_path, args.phylogenetic_raw_tree_path, args.tmp_path, args.seed,
                                   args.tree_search_software, args.outgroup, args.num_of_bootstrap_iterations, args.cpu)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
