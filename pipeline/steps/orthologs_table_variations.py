import itertools
import argparse
from pathlib import Path
import sys
import pandas as pd
import re
from ete3 import orthoxml
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
from pipeline.auxiliaries.general_utils import str_to_bool
from pipeline.auxiliaries import consts

ORPHANS_FILENAME_GENOME_NAME_PATTERN = re.compile('(.+)_orphans.txt')


def parse_species_name(logger, name, qfo_benchmark=False):
    if not qfo_benchmark:
        return name, None, None

    try:
        database_name, ncbi_tax_id, _ = name.split('_')
        return None, ncbi_tax_id, database_name
    except Exception:
        logger.exception(f"Failed parsing species name {name}, using full name")
        return name, None, None


def parse_gene_name(logger, name, qfo_benchmark=False):
    if not qfo_benchmark:
        return name

    try:
        prot_id = name.split('|')[1]
        return prot_id
    except Exception:
        logger.exception(f"Failed parsing gene name {name}, using full name")
        return name


def fix_orthoxml_output_file(orthoxml_file_path):
    with open(orthoxml_file_path, 'r') as oxml_file:
        content = oxml_file.read()

    # remove b'' prefixes
    pattern = re.compile(r'b\'\"([\w.:|-]+)\"\'')
    fixed_content = pattern.sub(r'"\1"', content)

    # remove 'ortho:' prefixes
    fixed_content = fixed_content.replace('ortho:', '')

    # add xmlns attribute to root xml tag
    fixed_content = fixed_content.replace('<orthoXML ', r'<orthoXML xmlns="http://orthoXML.org/2011/" ')

    with open(orthoxml_file_path, 'w') as oxml_file:
        oxml_file.write(fixed_content)


def build_orthoxml_and_tsv_output(logger, all_clusters_df, output_dir, qfo_benchmark=False):
    gene_name_to_gene_id = {}
    gene_name_to_gene_prot_id = {}

    # Creates an empty orthoXML object
    oxml = orthoxml.orthoXML(version="0.3", origin="Microbializer", originVersion="1")

    # Add the species (and their genes) to the orthoXML document
    next_id = 1
    for strain_column in all_clusters_df.columns[1:]:
        species_name, ncbi_tax_id, database_name = parse_species_name(logger, strain_column, qfo_benchmark)
        species_xml = orthoxml.species(name=species_name, NCBITaxId=ncbi_tax_id)
        database_xml = orthoxml.database(name=database_name)
        genes_xml = orthoxml.genes()

        all_strain_genes = all_clusters_df[strain_column]
        for og_index, og_strain_genes in all_strain_genes.items():
            if pd.isna(og_strain_genes):
                continue

            try:
                for gene_name in og_strain_genes.split(';'):
                    prot_id = parse_gene_name(logger, gene_name, qfo_benchmark)
                    gene_xml = orthoxml.gene(id=str(next_id), protId=prot_id)
                    genes_xml.add_gene(gene_xml)

                    gene_name_to_gene_id[gene_name] = next_id
                    gene_name_to_gene_prot_id[gene_name] = prot_id

                    next_id += 1
            except Exception as e:
                raise BrokenPipeError(f"Failed parsing gene names for column {strain_column} and OG {og_index}: {e}")

        database_xml.set_genes(genes_xml)
        species_xml.add_database(database_xml)
        oxml.add_species(species_xml)

    # Add an ortho group container to the orthoXML document
    groups_xml = orthoxml.groups()
    oxml.set_groups(groups_xml)

    # Add ortholog groups to the orthoXML document + construct a list of all ortholog groups
    ortholog_groups = []
    for index, group_row in all_clusters_df.iterrows():
        og_xml = orthoxml.group(id=group_row['OG_name'])
        og_list = []
        for strain_genes in group_row[1:]:
            if pd.isna(strain_genes):
                continue

            strain_genes = strain_genes.split(';')
            if len(strain_genes) == 1:
                gene = strain_genes[0]
                gene_ref_xml = orthoxml.geneRef(gene_name_to_gene_id[gene])
                og_xml.add_geneRef(gene_ref_xml)
                og_list.append([gene_name_to_gene_prot_id[gene]])
            else:
                paralog_group_xml = orthoxml.group()
                paralog_group_list = []
                for gene in strain_genes:
                    gene_ref_xml = orthoxml.geneRef(gene_name_to_gene_id[gene])
                    paralog_group_xml.add_geneRef(gene_ref_xml)
                    paralog_group_list.append(gene_name_to_gene_prot_id[gene])
                og_xml.add_paralogGroup(paralog_group_xml)
                og_list.append(paralog_group_list)

        groups_xml.add_orthologGroup(og_xml)

        if consts.OUTPUT_TSV_OF_ORTHOLOGS_PAIRS:
            ortholog_groups.append(og_list)

    # export orthoXML document to output_file
    orthoxml_output_file_path = output_dir / 'orthogroups.orthoxml'
    with open(orthoxml_output_file_path, "w") as oxml_file:
        oxml.export(oxml_file, level=0)

    fix_orthoxml_output_file(orthoxml_output_file_path)
    logger.info(f'Created orthoxml output at {orthoxml_output_file_path}')

    # write to tsv output all ortholog pairs
    if consts.OUTPUT_TSV_OF_ORTHOLOGS_PAIRS:
        tsv_output_file_path = output_dir / 'ortholog_pairs.tsv'
        with open(tsv_output_file_path, "w") as tsv_file:
            for og_list in ortholog_groups:
                for strain1_genes, strain2_genes in itertools.combinations(og_list, 2):
                    for strain1_gene, strain2_gene in itertools.product(strain1_genes, strain2_genes):
                        tsv_file.write(f'{strain1_gene}\t{strain2_gene}\n')
        logger.info(f'Created tsv orthologs pairs output at {tsv_output_file_path}')


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
        # Hierarchical clustering linkage (you can try different metrics/methods)
        linkage_matrix = linkage(pdist(binary_df, metric='hamming'), method='average')

        # Plot with seaborn's clustermap (just cluster rows, not columns)
        g = sns.clustermap(
            binary_df,
            row_linkage=linkage_matrix,
            col_cluster=False,  # disable clustering on orthogroups
            cmap="Blues",  # black = presence, white = absence (or reverse)
            figsize=compute_figsize(*binary_df.shape),
            xticklabels=False,
            yticklabels=True
        )

        # Apply smaller font size to y-axis labels
        g.ax_heatmap.tick_params(axis='y', labelsize=compute_fontsize(binary_df.shape[0]))

        # Remove the x-axis label/title and color legend
        g.ax_heatmap.set_xlabel("")  # remove x-axis label
        g.ax_heatmap.set_title("")  # just in case a title sneaks in
        g.cax.set_visible(False)

        g.fig.subplots_adjust(top=0.90)

        plt.tight_layout()
        g.savefig(map_png_path, dpi=600, bbox_inches='tight')
        g.savefig(map_svg_path, dpi=600, bbox_inches='tight')
        plt.close()

        logger.info(f'Created presence/absence map at {map_png_path} and {map_svg_path}')
    except Exception as e:
        logger.exception(f"Error creating presence/absence map: {e}")


def cluster_strains_by_orthogroups(logger, binary_df, output_dir):
    n_strains = binary_df.shape[0]

    if n_strains < 30:
        logger.info("Not enough strains to run umap+strains clustering, skipping.")
        return

    def compute_figsize(n, base=10, scale=0.015, max_size=30):
        side = min(base + n * scale, max_size)
        return side, side

    strains_cluster_csv_path = output_dir / 'strain_cluster_mapping.csv'
    strains_clusters_png_path = output_dir / 'strain_clusters_by_orthogroups.png'
    strains_clusters_svg_path = output_dir / 'strain_clusters_by_orthogroups.svg'

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

        plt.title("Clustering of Strains by Orthogroups presence/absence pattern", fontsize=18)
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


def create_orthogroups_variations(logger, orthologs_table_path, output_dir, qfo_benchmark):
    orthogroups_df = pd.read_csv(orthologs_table_path, dtype=str)
    create_phyletic_pattern(logger, orthogroups_df, output_dir)

    # Transpose to make rows = strains, columns = orthogroups
    binary_df = orthogroups_df.set_index('OG_name').notna().astype(int).T
    plot_presence_absence_matrix(logger, binary_df, output_dir)

    strain_clusters_dir = output_dir / 'strain_clusters'
    strain_clusters_dir.mkdir(parents=True, exist_ok=True)
    cluster_strains_by_orthogroups(logger, binary_df, strain_clusters_dir)

    build_orthoxml_and_tsv_output(logger, orthogroups_df, output_dir, qfo_benchmark)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthologs_table_path', type=Path, help='path to the orthologs table (input)')
    parser.add_argument('output_dir', type=Path, help='path to the output dir')
    parser.add_argument('--qfo_benchmark', help='whether the output OrthoXml should be in QfO benchmark format',
                        type=str_to_bool)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, create_orthogroups_variations, args.orthologs_table_path, args.output_dir, args.qfo_benchmark)
