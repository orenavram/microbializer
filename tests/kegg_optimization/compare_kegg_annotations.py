from pathlib import Path
import pandas as pd


BASE_PATH = Path(r"C:\repos\microbializer\tests\kegg_optimization\salmonella_300_c_001")
USE_ALL_GENES_PATH = BASE_PATH / "use_all_genes_of_og" / "og_kegg.csv"
USE_FIRST_GENE_PATH = BASE_PATH / "use_first_gene_of_og" / "og_kegg.csv"
USE_CONSENSUS_PATH = BASE_PATH / "use_consensus_of_og" / "og_kegg.csv"


def parse_knums(knum):
    if pd.isna(knum):
        return set()
    return set(knum.split(';'))


def compute_metrics(row):
    true_labels = row['knum_parsed_true']
    approx_labels = row['knum_parsed_approx']

    tp = len(true_labels & approx_labels)  # Intersection (correct predictions)
    fp = len(approx_labels - true_labels)  # Predicted but not in true
    fn = len(true_labels - approx_labels)  # Missing true labels

    return pd.Series({'TP': tp, 'FP': fp, 'FN': fn})


def compare(true_df, optimized_df, optimization_mode):
    true_df['knum_parsed'] = true_df['knum'].apply(parse_knums)
    optimized_df['knum_parsed'] = optimized_df['knum'].apply(parse_knums)

    # Merge DataFrames on OG_name
    merged_df = pd.merge(true_df, optimized_df, on='OG_name', suffixes=('_true', '_approx'))

    metrics_df = merged_df.apply(compute_metrics, axis=1)

    # Compute overall precision, recall, and F1-score
    TP_total = metrics_df['TP'].sum()
    FP_total = metrics_df['FP'].sum()
    FN_total = metrics_df['FN'].sum()

    precision = TP_total / (TP_total + FP_total) if (TP_total + FP_total) > 0 else 0
    recall = TP_total / (TP_total + FN_total) if (TP_total + FN_total) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    metrics = {
        'optimization_mode': optimization_mode,
        'precision': round(precision, 2),
        'recall': round(recall, 2),
        'f1_score': round(f1_score, 2)
    }

    metrics_df = pd.DataFrame(metrics, index=[0])
    return metrics_df


def main():
    all_genes_df = pd.read_csv(USE_ALL_GENES_PATH)
    first_gene_df = pd.read_csv(USE_FIRST_GENE_PATH)
    consensus_df = pd.read_csv(USE_CONSENSUS_PATH)

    metrics_all = compare(all_genes_df, all_genes_df, 'all_genes_of_og')
    metrics_first = compare(all_genes_df, first_gene_df, 'first_gene_of_og')
    metrics_consensus = compare(all_genes_df, consensus_df, 'consensus_of_og')

    metrics = pd.concat([metrics_all, metrics_first, metrics_consensus], axis=0)
    metrics.to_csv(BASE_PATH / 'metrics.csv', index=False)


if __name__ == "__main__":
    main()
