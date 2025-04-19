import argparse
from pathlib import Path
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
import numpy as np
import umap
from sklearn.cluster import HDBSCAN

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.logic_utils import count_strains_and_genes_in_ogs


def create_phyletic_pattern(logger, orthogroups_df, output_dir):
    phyletic_patterns_str = ''
    strain_names = list(orthogroups_df.columns[1:])
    for strain_name in strain_names:
        phyletic_pattern = ''.join(pd.notnull(orthogroups_df[strain_name]).astype(int).astype(str))
        phyletic_patterns_str += f'>{strain_name}\n{phyletic_pattern}\n'

    phyletic_patterns_path = output_dir / 'phyletic_pattern.fas'
    with open(phyletic_patterns_path, 'w') as f:
        f.write(phyletic_patterns_str)
    logger.info(f'Created phyletic pattern at {phyletic_patterns_path}')


def plot_presence_absence_matrix(logger, binary_df, output_dir):
    def compute_figsize(n_rows, n_cols, max_width=30, max_height=20, min_width=8, min_height=6):
        # Logarithmic scaling with caps
        width = min(max(min_width, np.log10(n_cols) * 5), max_width)
        height = min(max(min_height, np.log10(n_rows) * 2.5), max_height)
        return width, height

    def compute_fontsize(n_rows, min_font=1, max_font=10):
        return max(min_font, min(max_font, 200 / n_rows))  # smaller than before

    map_png_path = output_dir / 'phyletic_pattern.png'
    map_svg_path = output_dir / 'phyletic_pattern.svg'

    try:
        n_genomes = binary_df.shape[0]
        n_orthogroups = binary_df.shape[1]

        if n_genomes < 3:
            logger.info("Not enough strains to run clustering. Will produce heatmap instead.")
            linkage_matrix = None
            row_cluster = False
        else:
            # Hierarchical clustering linkage
            linkage_matrix = linkage(pdist(binary_df, metric='hamming'), method='average')
            row_cluster = True

        # Plot with seaborn's clustermap (just cluster rows, not columns)
        g = sns.clustermap(
            binary_df,
            row_linkage=linkage_matrix,
            row_cluster=row_cluster,
            col_cluster=False,  # disable clustering on orthogroups
            cmap="Blues",  # black = presence, white = absence (or reverse)
            figsize=compute_figsize(n_genomes, n_orthogroups),
        )

        # Remove the x-axis label, title and color legend
        g.ax_heatmap.set_xlabel("")  # remove x-axis label
        g.ax_heatmap.set_title("")  # just in case a title sneaks in
        g.cax.set_visible(False)

        plt.tight_layout()
        g.savefig(map_png_path, dpi=600, bbox_inches='tight')
        g.savefig(map_svg_path, dpi=600, bbox_inches='tight')
        plt.close()

        logger.info(f'Created presence/absence map at {map_png_path} and {map_svg_path}')
    except Exception as e:
        logger.exception(f"Error creating presence/absence map: {e}")


def cluster_strains_by_orthogroups(logger, binary_df, output_dir):
    def compute_figsize(n, base=10, scale=0.015, max_size=30):
        side = min(base + n * scale, max_size)
        return side, side

    strains_cluster_csv_path = output_dir / 'strain_cluster_mapping.csv'
    strains_clusters_png_path = output_dir / 'strain_clusters_by_orthogroups.png'
    strains_clusters_svg_path = output_dir / 'strain_clusters_by_orthogroups.svg'

    n_strains = binary_df.shape[0]
    if n_strains < 30:
        logger.info("Not enough strains to run umap and strains clustering, skipping.")
        open(strains_cluster_csv_path, 'w').close()
        open(strains_clusters_png_path, 'w').close()
        return

    try:
        # Step 1: UMAP dimensionality reduction (good with Hamming distance for binary data)
        reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='hamming', random_state=42, n_jobs=1)
        embedding = reducer.fit_transform(binary_df)

        # Step 2: HDBSCAN clustering (no need to specify number of clusters)
        min_cluster_size = max(5, int(n_strains * 0.01))  # Dynamically set min_cluster_size
        clusterer = HDBSCAN(min_cluster_size=min_cluster_size, metric='hamming')
        labels = clusterer.fit_predict(binary_df)

        # 3. Save strain-cluster mapping to CSV
        strain_cluster_df = pd.DataFrame({
            'strain': binary_df.index,
            'cluster': labels
        })
        strain_cluster_df.to_csv(strains_cluster_csv_path, index=False)

        # 4. Build color palette with unique color for -1
        unique_labels = np.unique(labels)
        n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
        base_palette = sns.color_palette("husl", n_colors=max(10, n_clusters))

        # Define palette: normal clusters get colors, noise gets gray
        palette = {}
        color_idx = 0
        for label in unique_labels:
            if label == -1:
                palette[label] = (0.6, 0.6, 0.6)  # gray for noise
            else:
                palette[label] = base_palette[color_idx]
                color_idx += 1

        # Step 5: Plot
        plt.figure(figsize=compute_figsize(n_strains))
        sns.scatterplot(
            x=embedding[:, 0],
            y=embedding[:, 1],
            hue=labels,
            palette=palette,
            s=50,
            legend='full'
        )

        plt.title("Clustering of Strains by Orthogroups presence/absence pattern", fontsize=20)
        plt.xlabel("UMAP-1")
        plt.ylabel("UMAP-2")
        plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        # Save
        plt.savefig(strains_clusters_png_path, dpi=600)
        plt.savefig(strains_clusters_svg_path)
        plt.close()

        logger.info(f'Created strains clusters by orthogroups at {strains_clusters_png_path} and {strains_clusters_svg_path}')
    except Exception as e:
        logger.exception(f"Error creating strains clusters by orthogroups: {e}")


def count_and_plot_orthogroups_sizes(logger, final_orthogroups_file_path, output_dir):
    final_orthologs_table_df = pd.read_csv(final_orthogroups_file_path, dtype=str)
    orthogroups_sizes_df = count_strains_and_genes_in_ogs(final_orthologs_table_df)
    orthogroups_sizes_df = orthogroups_sizes_df.rename(columns={'strains_count': 'OG size (number of genomes)',
                                                                'genes_count': 'OG size (total number of genes)'})
    orthogroups_sizes_df.to_csv(output_dir / 'orthogroups_sizes.csv', index=False)
    logger.info(f'Created orthogroups sizes table at {output_dir / "orthogroups_sizes.csv"}')

    group_sizes = orthogroups_sizes_df.set_index('OG_name')['OG size (number of genomes)']
    sns.histplot(x=group_sizes, discrete=True)
    if len(np.unique(group_sizes)) < 10:
        plt.xticks(np.unique(group_sizes))
    plt.title('Orthogroups sizes distribution', fontsize=25, loc='center', wrap=True)
    plt.xlabel('OG size (number of genomes)', fontsize=20)
    plt.ylabel('Count of OGs of each OG size', fontsize=20)
    plt.tight_layout()
    plt.savefig(output_dir / 'orthogroups_sizes.png', dpi=600)
    plt.savefig(output_dir / 'orthogroups_sizes.svg', dpi=600)
    plt.clf()
    logger.info(f'Created orthogroups sizes histogram at {output_dir / "orthogroups_sizes.png"}')


def create_simplified_orthogroups_table_for_results_page(logger, orthogroups_df, output_dir):
    orthogroups_df = orthogroups_df.set_index('OG_name')

    # Start with a new DataFrame with the same shape as the original, filled with 0.
    converted_df = pd.DataFrame(0, index=orthogroups_df.index, columns=orthogroups_df.columns)

    not_nan_mask = orthogroups_df.notna()
    semicolon_mask = orthogroups_df.apply(lambda col: col.str.contains(';', na=False))

    converted_df[not_nan_mask] = 1
    converted_df[semicolon_mask] = 2

    results_page_orthogroups_path = output_dir / 'orthogroups_results_page.csv'
    converted_df.to_csv(results_page_orthogroups_path)
    logger.info(f'Created simplified orthogroups table for results page at {results_page_orthogroups_path}')


def create_orthogroups_visualizations(logger, orthologs_table_path, output_dir, tmp_dir):
    count_and_plot_orthogroups_sizes(logger, orthologs_table_path, output_dir)

    orthogroups_df = pd.read_csv(orthologs_table_path, dtype=str)
    create_simplified_orthogroups_table_for_results_page(logger, orthogroups_df, tmp_dir)
    create_phyletic_pattern(logger, orthogroups_df, output_dir)

    # Transpose to make rows = strains, columns = orthogroups
    binary_df = orthogroups_df.set_index('OG_name').notna().astype(int).T

    plot_presence_absence_matrix(logger, binary_df, output_dir)
    cluster_strains_by_orthogroups(logger, binary_df, output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', type=Path, help='path to the orthologs table (input)')
    parser.add_argument('output_dir', type=Path, help='path to the output dir')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, create_orthogroups_visualizations, args.orthologs_table_path, args.output_dir, args.logs_dir)
