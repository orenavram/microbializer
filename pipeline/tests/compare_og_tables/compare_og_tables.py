import os
import pandas as pd
from sklearn import metrics
from collections import Counter

BASE_PATH = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\73_ecoli"
OG_TABLE_NO_OPTIMIZE = os.path.join(BASE_PATH, "no_optimize", "final_orthologs_table.csv")
OG_TABLE_OPTIMIZE_5 = os.path.join(BASE_PATH, "optimize_5_clusters", "final_orthologs_table.csv")
OG_TABLE_OPTIMIZE_10 = os.path.join(BASE_PATH, "optimize_10_clusters", "final_orthologs_table.csv")


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
    true_labels, true_genes = convert_og_table_to_labels(OG_TABLE_NO_OPTIMIZE)
    pred_labels_5, pred_genes_5 = convert_og_table_to_labels(OG_TABLE_OPTIMIZE_5)
    pred_labels_10, pred_genes_10 = convert_og_table_to_labels(OG_TABLE_OPTIMIZE_10)

    if len(true_labels) != len(pred_labels_5):
        genes_in_true_not_in_5 = set(true_genes).difference(set(pred_genes_5))
        print(f"Genes in true but not in pred_5: {genes_in_true_not_in_5}")
        genes_in_5_not_in_true = set(pred_genes_5).difference(set(true_genes))
        print(f"Genes in pred_5 but not in true: {genes_in_5_not_in_true}")
        true_duplicates = [item for item, count in Counter(true_genes).items() if count > 1]
        pred_5_duplicates = [item for item, count in Counter(pred_genes_5).items() if count > 1]
        print(f"Duplicated in true: {true_duplicates}, duplicated in pred_5: {pred_5_duplicates}")
        raise ValueError("The number of genes in the ortholog tables true and pred_5 is not equal")
    if len(true_labels) != len(pred_labels_10):
        pred_10_duplicated = [item for item, count in Counter(pred_genes_10).items() if count > 1]
        print(f"Duplicated in pred_10: {pred_10_duplicated}")
        raise ValueError("The number of genes in the ortholog tables true and pred_10 is not equal")

    compare_clusterings(true_labels, pred_labels_5, os.path.join(BASE_PATH, "comparison_scores_true_vs_5.csv"))
    compare_clusterings(true_labels, pred_labels_10, os.path.join(BASE_PATH, "comparison_scores_true_vs_10.csv"))
    compare_clusterings(pred_labels_5, pred_labels_10, os.path.join(BASE_PATH, "comparison_scores_5_vs_10.csv"))


if __name__ == '__main__':
    main()
