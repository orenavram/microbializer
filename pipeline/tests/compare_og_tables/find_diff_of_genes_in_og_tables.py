import pandas as pd
from collections import Counter

PATH_1 = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\4_genomes\optimize_1_mmseqs_command\putative_orthologs_table.txt"
PATH_2 = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\4_genomes\no_optimize\putative_orthologs_table.txt"

df_1 = pd.read_csv(PATH_1).drop(columns=['OG_name'])
all_values_1 = df_1.to_numpy().ravel()
all_values_1 = all_values_1[~(pd.isna(all_values_1) | (all_values_1 == ""))]
all_genes_1 = [item for sublist in all_values_1 for item in sublist.split(';')]

df_2 = pd.read_csv(PATH_2).drop(columns=['OG_name'])
all_values_2 = df_2.to_numpy().ravel()
all_values_2 = all_values_2[~(pd.isna(all_values_2) | (all_values_2 == ""))]
all_genes_2 = [item for sublist in all_values_2 for item in sublist.split(';')]

print(f"Number of genes in PATH_1: {len(all_genes_1)}, Number of genes in PATH_2: {len(all_genes_2)}")

diff_1_to_2 = [gene for gene in all_genes_1 if gene not in all_genes_2]
diff_2_to_1 = [gene for gene in all_genes_2 if gene not in all_genes_1]

print(f"Number of genes in PATH_1 but not in PATH_2: {len(diff_1_to_2)}, genes: {diff_1_to_2}")
print(f"Number of genes in PATH_2 but not in PATH_1: {len(diff_2_to_1)}, genes: {diff_2_to_1}")
