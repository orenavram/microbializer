import os
import pandas as pd
from sklearn import metrics

BASE_PATH = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\example_data_10_genomes"
OG_TABLE_NO_OPTIMIZE = os.path.join(BASE_PATH, "final_orthologs_table_no_optimize.csv")
OG_TABLE_OPTIMIZE = os.path.join(BASE_PATH, "final_orthologs_table_optimize.csv")
OUTPUT_PATH = os.path.join(BASE_PATH, "comparison_scores.csv")


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

    return labels_list


def main():
    true_labels = convert_og_table_to_labels(OG_TABLE_NO_OPTIMIZE)
    pred_labels = convert_og_table_to_labels(OG_TABLE_OPTIMIZE)

    if len(true_labels) != len(pred_labels):
        raise ValueError("The number of genes in the two ortholog tables is not equal")

    compare_clusterings(true_labels, pred_labels, OUTPUT_PATH)


if __name__ == '__main__':
    main()
