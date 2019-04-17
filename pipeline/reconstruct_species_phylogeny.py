"""

"""

# python /bioseq/microbializer/pipeline/reconstruct_species_phylogeny.py /bioseq/data/results/microbializer/155542823177458857633443576357/outputs/15_aligned_core_proteome/aligned_core_proteome.fas /bioseq/data/results/microbializer/155542823177458857633443576357/outputs/16_species_phylogeny/species_tree.txt --model PROTGAMMAILG --num_of_bootstrap_iterations 100 --cpu 3 --num_of_bootstrap_iterations 5 --root


def generate_phylogenetic_tree(msa_path, phylogenetic_tree_path, seed, model, num_of_bootstrap_iterations, root, num_of_cpus):
    import os
    wd, final_tree_name = os.path.split(phylogenetic_tree_path)

    intermediate_tree_path = phylogenetic_tree_path
    phylogenetic_tree_name = 'unrooted_species_tree.txt'
    rooted_phylogenetic_tree_name = 'rooted_species_tree.txt'

    import subprocess
    #e.g.: raxmlHPC-PTHREADS-SSE3 -m PROTGAMMAILG -p 12345 -s /bioseq/data/results/microbializer/155542823177458857633443576357/outputs/15_aligned_core_proteome/aligned_core_proteome.fas -n unrooted_species_tree.txt -w /bioseq/data/results/microbializer/155542823177458857633443576357/outputs/16_species_phylogeny -T 3 -f a -x 12345 -N 5
    cmd = f'raxmlHPC-PTHREADS-SSE3 -m {model} -p {seed} -s {msa_path} -n {phylogenetic_tree_name} -w {wd} -T {num_of_cpus}'
    if int(num_of_bootstrap_iterations)>0:
        cmd += ' ' + f'-f a -x {seed} -N {num_of_bootstrap_iterations}'
    logger.info(f'Reconstructing species phylogeny with RAxML. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)

    raxml_best_tree_output_path = os.path.join(wd, f'RAxML_bestTree.{phylogenetic_tree_name}')

    bootstrapped_raw_tree_path = os.path.join(wd, f'RAxML_bipartitions.{phylogenetic_tree_name}')

    if os.path.exists(bootstrapped_raw_tree_path):
        intermediate_tree_path = bootstrapped_raw_tree_path
    elif os.path.exists(raxml_best_tree_output_path):
        intermediate_tree_path = bootstrapped_raw_tree_path
    else:
        #TODO: make sure that the script of ploting the tree is aware to that there is no tree.
        logger.fatal(f'TREE WAS NOT GENERATED!!')
        pass

    if root and os.path.exists(intermediate_tree_path):
        # raxmlHPC-PTHREADS-SSE3 -m PROTGAMMAILG -n rooted_species_tree.txt -w /bioseq/data/results/microbializer/155542823177458857633443576357/outputs/16_species_phylogeny -t /bioseq/data/results/microbializer/155542823177458857633443576357/outputs/16_species_phylogeny/RAxML_bipartitions.unrooted_species_tree.txt -f I
        cmd = f'raxmlHPC-PTHREADS-SSE3  -m {model} -n {rooted_phylogenetic_tree_name} -w {wd} -t {intermediate_tree_path} -f I'
        logger.info(f'Rooting species tree. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)
        if os.path.exists(os.path.join(wd, f'RAxML_rootedTree.{rooted_phylogenetic_tree_name}')):
            intermediate_tree_path = os.path.join(wd, f'RAxML_rootedTree.{rooted_phylogenetic_tree_name}')
            logger.info(f'Rooting succeeded!')
        else:
            logger.fatal(f'Rooting failed! Returning the unrooted tree.')

    os.rename(intermediate_tree_path, phylogenetic_tree_path)

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
        parser.add_argument('phylogenetic_raw_tree_path', help='path to an output file in which the phylogenetic tree will be written')
        parser.add_argument('--seed',
                            help='RaxML seed parameter',
                            default=12345)
        parser.add_argument('--model',
                            help='RAxML reconstruction model parameter')
        parser.add_argument('--num_of_bootstrap_iterations', default=0)
        parser.add_argument('--root', action='store_true')
        parser.add_argument('--cpu', choices=[str(i) for i in range(1, 29)], default='1',
                            help='How many CPUs will be used? (for running in parallel mode). For further details see:\nhttps://support.nesi.org.nz/hc/en-gb/articles/115001854444-RAxML#parallel-versions')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        generate_phylogenetic_tree(args.msa_path, args.phylogenetic_raw_tree_path, args.seed,
                                   args.model, args.num_of_bootstrap_iterations, args.root, args.cpu)

