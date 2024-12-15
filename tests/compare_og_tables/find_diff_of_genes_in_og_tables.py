from pathlib import Path
import pandas as pd
from collections import Counter

BASE_PATH = Path(r"C:\temp")
PATH_1 = BASE_PATH / "putative_orthologs_table.txt"
PATH_2 = BASE_PATH / "orthogroups.csv"


df_1 = pd.read_csv(PATH_1).drop(columns=['OG_name'])
all_values_1 = df_1.to_numpy().ravel()
all_values_1 = all_values_1[~(pd.isna(all_values_1) | (all_values_1 == ""))]
all_genes_1 = set([item for sublist in all_values_1 for item in sublist.split(';')])

df_2 = pd.read_csv(PATH_2).drop(columns=['OG_name'])
all_values_2 = df_2.to_numpy().ravel()
all_values_2 = all_values_2[~(pd.isna(all_values_2) | (all_values_2 == ""))]
all_genes_2 = set([item for sublist in all_values_2 for item in sublist.split(';')])

print(f"Number of genes in PATH_1: {len(all_genes_1)}, Number of genes in PATH_2: {len(all_genes_2)}")

diff_1_to_2 = all_genes_1 - all_genes_2
with open(BASE_PATH / "diff_1_to_2.txt", 'w') as f:
    f.write("\n".join(diff_1_to_2))
diff_2_to_1 = all_genes_2 - all_genes_1
with open(BASE_PATH / "diff_2_to_1.txt", 'w') as f:
    f.write("\n".join(diff_2_to_1))

print(f"Number of genes in PATH_1 but not in PATH_2: {len(diff_1_to_2)}")
print(f"Number of genes in PATH_2 but not in PATH_1: {len(diff_2_to_1)}")
