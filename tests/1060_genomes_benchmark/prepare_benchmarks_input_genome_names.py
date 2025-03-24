import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

XANTHOMONAS_METADATA = SCRIPT_DIR / "ncbi_xanthomonas.tsv"
ECOLI_METADATA = SCRIPT_DIR / "ncbi_ecoli.tsv"

xanthomonas_df = pd.read_csv(XANTHOMONAS_METADATA, sep='\t')
ecoli_df = pd.read_csv(ECOLI_METADATA, sep='\t')

xanthomonas_accessions = xanthomonas_df['Assembly Accession'].tolist()
ecoli_accessions = ecoli_df['Assembly Accession'].tolist()

all_accessions = xanthomonas_accessions + ecoli_accessions
with open(SCRIPT_DIR / 'all_accessions.txt', 'w') as f:
    f.write('\n'.join(all_accessions))
