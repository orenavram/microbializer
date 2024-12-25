import pandas as pd

PATH_1 = r"C:\temp\use_dbs\Sflexneri_5_8401_vs_Ecoli_O104_H4_2009EL_2050.m8"
PATH_2 = r"C:\temp\original\Sflexneri_5_8401_vs_Ecoli_O104_H4_2009EL_2050.m8"


df_1 = pd.read_csv(PATH_1)
all_genes_1 = set(df_1['query'].tolist() + df_1['target'].tolist())

df_2 = pd.read_csv(PATH_2)
all_genes_2 = set(df_2['query'].tolist() + df_2['target'].tolist())

print(f"Number of genes in PATH_1: {len(all_genes_1)}, Number of genes in PATH_2: {len(all_genes_2)}")

diff_1_to_2 = all_genes_1 - all_genes_2
diff_2_to_1 = all_genes_2 - all_genes_1

print(f"Number of genes in PATH_1 but not in PATH_2: {len(diff_1_to_2)}, genes: {diff_1_to_2}")
print(f"Number of genes in PATH_2 but not in PATH_1: {len(diff_2_to_1)}, genes: {diff_2_to_1}")
