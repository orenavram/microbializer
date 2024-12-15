import os
import subprocess
import pandas as pd

DOWNLOADS_DIR = r"C:\repos\microbializer\benchmark_data\downloads"
XANTHOMONAS_METADATA = r"C:\repos\microbializer\benchmark_code\ncbi_xanthomonas.tsv"
ECOLI_METADATA = r"C:\repos\microbializer\benchmark_code\ncbi_ecoli.tsv"

xanthomonas_df = pd.read_csv(XANTHOMONAS_METADATA, sep='\t')
ecoli_df = pd.read_csv(ECOLI_METADATA, sep='\t')

xanthomonas_accessions = xanthomonas_df['Assembly Accession'].tolist()
ecoli_accessions = ecoli_df['Assembly Accession'].tolist()

for accession in xanthomonas_accessions + ecoli_accessions:
    download_path = os.path.join(DOWNLOADS_DIR, f"{accession}.zip")
    cmd = f'datasets download genome accession {accession} --include genome,protein,cds,gff3 --filename {download_path}'
    subprocess.run(cmd, shell=True)