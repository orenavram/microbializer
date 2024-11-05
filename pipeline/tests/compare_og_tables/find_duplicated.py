import pandas as pd
from collections import Counter

BASE_PATH = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\73_ecoli\optimize_10_clusters\final_orthologs_table.csv"

df = pd.read_csv(BASE_PATH).drop(columns=['OG_name'])
all_values = df.to_numpy().ravel()
all_values = all_values[~(pd.isna(all_values) | (all_values == ""))]
all_genes = [item for sublist in all_values for item in sublist.split(';')]

duplicated = [item for item, count in Counter(all_genes).items() if count > 1]
print(f"Number of duplicated genes: {len(duplicated)}, duplicated genes: {duplicated}")
