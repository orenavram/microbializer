import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from Bio import Phylo
import pylab as pb


def generate_boxplot(path_to_data, output_file_path, title='', xlabel='', ylabel='', dpi=300):
    data = np.loadtxt(path_to_data)
    fig = plt.figure(figsize=(10,10))

    ax = sns.violinplot(data, inner=None, color='lavender', cut=5)
    # ax = sns.boxplot(data, color='lavender')
    sns.swarmplot(data, color='dodgerblue')
    ax.set_xlabel(f'{xlabel}', fontdict={'fontsize': 20})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 20})
    ax.set_title(f'{title}', fontdict={'fontsize': 20})

    fig.savefig(output_file_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def generate_tree_plot(path_to_data, output_file_path):
    tree = Phylo.read(path_to_data, 'newick')
    Phylo.draw(tree)
    pb.savefig(output_file_path)


def generate_bar_plot(path_to_data, output_file_path, title='', xlabel='', ylabel='', dpi=300):
    data = np.loadtxt(path_to_data)

    fig = plt.figure(figsize=(10, 10))

    bins = np.arange(data.min(), data.max() + 2)

    ax = sns.distplot(data, bins=bins, kde=False, hist_kws={"rwidth":0.9,'edgecolor':'black', 'alpha':0.5})
    ax.set_xlabel(f'{xlabel}', fontdict={'fontsize': 20})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 20})
    ax.set_title(f'{title}', fontdict={'fontsize': 20})
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    fig.savefig(output_file_path, dpi=dpi, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    generate_bar_plot(
        '/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/13_groups_sizes_frequency/groups_sizes_frequency.txt',
        '/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/13_groups_sizes_frequency/groups_sizes_frequency.png',
        xlabel='\nOrthologs group size', ylabel='Counts\n', dpi=100)
    generate_boxplot('/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/14_orfs_statistics/orfs_gc_contents.txt',
                      '/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/14_orfs_statistics/orfs_gc_contents.png',
                     xlabel='\nGC content per genome', dpi=100)
    generate_boxplot('/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/14_orfs_statistics/orfs_counts.txt',
                      '/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/14_orfs_statistics/orfs_counts.png',
                     xlabel='\nORFs count per genome', dpi=100)
    generate_tree_plot('/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/17_species_phylogeny/species_tree.txt',
                      '/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/17_species_phylogeny/species_tree.png')
