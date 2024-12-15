import pandas as pd
from collections import Counter

BASE_PATH = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\4_genomes\optimize_1_mmseqs_command\putative_orthologs_table.txt"

df = pd.read_csv(BASE_PATH).drop(columns=['OG_name'])
all_values = df.to_numpy().ravel()
all_values = all_values[~(pd.isna(all_values) | (all_values == ""))]
all_genes = [item for sublist in all_values for item in sublist.split(';')]

duplicated = [item for item, count in Counter(all_genes).items() if count > 1]
print(f"Number of duplicated genes: {len(duplicated)} (out of total {len(all_genes)}), duplicated genes: {duplicated}")
