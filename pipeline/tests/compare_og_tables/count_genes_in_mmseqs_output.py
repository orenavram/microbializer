import pandas as pd

PATH_1 = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\4_genomes\no_optimize\05c_mmseqs_paralogs\NC_007606_vs_NC_007606.m8_filtered"
PATH_2 = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\4_genomes\optimize_1_mmseqs_command\05c_paralogs\NC_007606_vs_NC_007606_filtered.m8"

df_1 = pd.read_csv(PATH_1)
all_genes_1 = set(df_1['query'].tolist() + df_1['target'].tolist())

df_2 = pd.read_csv(PATH_2)
all_genes_2 = set(df_2['query'].tolist() + df_2['target'].tolist())

print(f"Number of genes in PATH_1: {len(all_genes_1)}, Number of genes in PATH_2: {len(all_genes_2)}")

diff_1_to_2 = all_genes_1 - all_genes_2
diff_2_to_1 = all_genes_2 - all_genes_1

print(f"Number of genes in PATH_1 but not in PATH_2: {len(diff_1_to_2)}, genes: {diff_1_to_2}")
print(f"Number of genes in PATH_2 but not in PATH_1: {len(diff_2_to_1)}, genes: {diff_2_to_1}")
