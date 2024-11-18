import pandas as pd

PATH_1 = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\4_genomes\no_optimize\putative_orthologs_table.txt"
PATH_2 = r"C:\repos\microbializer\pipeline\tests\compare_og_tables\4_genomes\optimize_1_mmseqs_command\putative_orthologs_table.txt"

df_1 = pd.read_csv(PATH_1)
all_genes_1 = set(df_1['query'].tolist() + df_1['target'].tolist())

df_2 = pd.read_csv(PATH_2)
all_genes_2 = set(df_2['query'].tolist() + df_2['target'].tolist())

print(f"Number of genes in PATH_1: {len(all_genes_1)}, Number of genes in PATH_2: {len(all_genes_2)}")
