import shutil
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

DOWNLOADS_DIR = SCRIPT_DIR / 'downloads'
EXTRACTED_DIR = SCRIPT_DIR / 'extracted'
ORFS_DIR = SCRIPT_DIR / 'orfs'
NCBI_METADATA = SCRIPT_DIR / 'ncbi_dataset_only_organism.tsv'

df = pd.read_csv(NCBI_METADATA)
accession_to_organism = dict(zip(df['Assembly Accession'], df['Organism Name']))

for folder in EXTRACTED_DIR.iterdir():
    accession = folder.name
    # Construct full folder path
    orfs_file_path = folder / 'ncbi_dataset' / 'data' / accession / 'cds_from_genomic.fna'

    # Copy orfs file to ORFS_DIR
    shutil.copy(orfs_file_path, ORFS_DIR / f'{accession_to_organism[accession]}.fna')
