import matplotlib

matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
import pylab as pb


def generate_violinplot(path_to_data, output_file_path, title='', xlabel='', ylabel='', dpi=300):
    data = np.atleast_1d(np.loadtxt(path_to_data))
    fig = plt.figure(figsize=(10, 10))

    ax = sns.violinplot(x=data, inner=None, color='lavender', cut=5)
    # ax = sns.boxplot(data, color='lavender')
    sns.swarmplot(x=data, color='dodgerblue')
    # ax.set_xlim(left=0)
    ax.set_xlabel(f'{xlabel}', fontdict={'fontsize': 20})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 20})
    ax.set_title(f'{title}', fontdict={'fontsize': 20})

    fig.savefig(output_file_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def generate_tree_plot(path_to_data, output_file_path):
    tree = Phylo.read(path_to_data, 'newick')
    Phylo.draw(tree)  # TODO: make figure proportional to number of species
    pb.savefig(output_file_path, bbox_inches='tight', dpi=300)


def generate_bar_plot(path_to_data, output_file_path, xlabel='', ylabel='', dpi=300):
    data = np.loadtxt(path_to_data, dtype=int)

    fig = plt.figure(figsize=(10, 10))

    # bins = np.arange(data.min(), data.max() + 2)

    ax = sns.countplot(x=data, color="C0")
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    visible_bins = max(1, len(np.unique(data)) // 20)
    print(f'len(np.unique(data))={len(np.unique(data))}')
    for ind, label in enumerate(ax.get_xticklabels()):
        if ind % visible_bins == 0:  # every $visible_bins bins, label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)
    ax.set_xlabel(f'{xlabel}', fontdict={'fontsize': 20})
    ax.set_ylabel(f'{ylabel}', fontdict={'fontsize': 20})

    # ax.set_title(f'{title}', fontdict={'fontsize': 20})
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    fig.savefig(output_file_path, dpi=dpi, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    dataset_name = '29 GammaProteoBacteria (no Philadelphia) dataset'

    generate_bar_plot(
        f'/Users/Oren/Dropbox/Projects/Sweeps Project/microbilizer_web/data/{dataset_name}/19_groups_sizes_frequency/groups_sizes_frequency.txt',
        f'/Users/Oren/Dropbox/Projects/microbializerPower/html/pics/{dataset_name.split()[0]}_groups_sizes_frequency.png',
        xlabel='Orthologs group size', ylabel='Counts')
    generate_violinplot(
        f'/Users/Oren/Dropbox/Projects/Sweeps Project/microbilizer_web/data/{dataset_name}/20_orfs_plots/orfs_gc_contents.txt',
        f'/Users/Oren/Dropbox/Projects/microbializerPower/html/pics/{dataset_name.split()[0]}_orfs_gc_contents.png',
        xlabel='GC content per genome')
    generate_violinplot(
        f'/Users/Oren/Dropbox/Projects/Sweeps Project/microbilizer_web/data/{dataset_name}/20_orfs_plots/orfs_counts.txt',
        f'/Users/Oren/Dropbox/Projects/microbializerPower/html/pics/{dataset_name.split()[0]}_orfs_counts.png',
        xlabel='ORFs count per genome')
    # generate_tree_plot('/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/17_species_phylogeny/final_species_tree.txt',
    #                  '/Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/17_species_phylogeny/final_species_tree.png')
