import os
import json
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def plot_genomes_histogram(data, output_dir, output_file_name, title, xlabel):
    # data is expected to be: {'genome1': 54, 'genome2': 20, ...}

    with open(os.path.join(output_dir, f'{output_file_name}.json'), 'w') as f:
        json.dump(data, f)

    output_df = pd.DataFrame.from_dict(data, orient='index', columns=[title])
    output_df.index.name = 'Genome'
    output_df.to_csv(os.path.join(output_dir, f'{output_file_name}.csv'))

    histogram = sns.histplot(output_df, x=title, kde=True)
    plt.title(f'Distribution of {title}')
    plt.xlabel(xlabel)
    plt.ylabel('Genomes count')
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig(os.path.join(output_dir, f'{output_file_name}.png'))

    plt.clf()


def main():
    with open(os.path.join(BASE_DIR, 'orfs_gc_contents.json'), 'r') as f:
        data = json.load(f)

    plot_genomes_histogram(data, BASE_DIR, 'gc_content', 'ORFs GC content', 'GC content per genome')


if __name__ == '__main__':
    main()
