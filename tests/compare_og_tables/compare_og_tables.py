import os
import pandas as pd
from sklearn import metrics

BASE_PATH = r"C:\Users\TalPNB22\Downloads"
OG_TABLE_NO_CLUSTER = os.path.join(BASE_PATH, "73", "orthogroups.csv")
OG_TABLE_CLUSTER = os.path.join(BASE_PATH, "10", "orthogroups.csv")


def compare_clusterings(true_labels, pred_labels, output_path):
    # https://scikit-learn.org/stable/modules/clustering.html#rand-index
    adjusted_rand_score = metrics.adjusted_rand_score(true_labels, pred_labels)

    # https://scikit-learn.org/stable/modules/clustering.html#mutual-information-based-scores
    adjusted_mutual_info_score = metrics.adjusted_mutual_info_score(true_labels, pred_labels)

    # https://scikit-learn.org/stable/modules/clustering.html#homogeneity-completeness-and-v-measure
    homogeneity, completeness, v_measure = metrics.homogeneity_completeness_v_measure(true_labels, pred_labels)

    # https://scikit-learn.org/stable/modules/clustering.html#fowlkes-mallows-scores
    fowlkes_mallows_score = metrics.fowlkes_mallows_score(true_labels, pred_labels)

    scores = {
        'adjusted_rand_score': adjusted_rand_score,
        'adjusted_mutual_info_score': adjusted_mutual_info_score,
        'homogeneity': homogeneity,
        'completeness': completeness,
        'v_measure': v_measure,
        'fowlkes_mallows_score': fowlkes_mallows_score
    }

    scores_df = pd.DataFrame([scores])
    scores_df.to_csv(output_path, index=False)


def convert_og_table_to_labels(og_table_path):
    df = pd.read_csv(og_table_path)

    # 1. Melt the DataFrame to get all genome entries in a single column
    melted_df = df.melt(id_vars='OG_name', value_name='genes', var_name='genome').dropna()

    # 2. Split the 'genes' column by ';' and explode into individual rows
    exploded_df = melted_df.assign(genes=melted_df['genes'].str.split(';')).explode('genes')

    # 3. Sort the DataFrame by gene names
    sorted_df = exploded_df[['genes', 'OG_name']].sort_values(by='genes')

    # 4. Extract the labels list sorted by gene names
    labels_list = sorted_df['OG_name'].tolist()

    genes_list = sorted_df['genes'].tolist()

    return labels_list, genes_list


def main():
    no_cluster_labels, no_cluster_genes = convert_og_table_to_labels(OG_TABLE_NO_CLUSTER)
    cluster_labels, cluster_genes = convert_og_table_to_labels(OG_TABLE_CLUSTER)

    if len(no_cluster_labels) != len(cluster_labels):
        genes_difference = set(no_cluster_genes).difference(set(cluster_genes))
        genes_difference.add(set(cluster_genes).difference(set(no_cluster_genes)))
        print(f"Genes in no_optimize but not in optimize or vice versa: {genes_difference}")
        return

    compare_clusterings(no_cluster_labels, cluster_labels, os.path.join(BASE_PATH, "comparison_scores.csv"))


if __name__ == '__main__':
    main()