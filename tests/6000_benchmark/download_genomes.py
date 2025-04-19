import os
import subprocess
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

DOWNLOADS_DIR = SCRIPT_DIR / 'downloads'
DOWNLOADS_DIR.mkdir(parents=True, exist_ok=True)

# This file was downloaded from NCBI Genome database on 19.4.2025 using the following filters:
# Search for all bacteria genomes
# include only: reference genomes, annotated genomes by NCBI RefSeq
# exclude: atypical, MAGs, multi-isolate projects
METADATA = SCRIPT_DIR / 'ncbi_dataset.tsv'

df = pd.read_csv(METADATA, sep='\t')
accessions = df['Assembly Accession'].tolist()

for accession in accessions:
    download_path = os.path.join(DOWNLOADS_DIR, f"{accession}.zip")
    cmd = f'datasets download genome accession {accession} --include cds --filename {download_path}'
    subprocess.run(cmd, shell=True)
