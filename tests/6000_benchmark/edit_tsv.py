from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
NCBI_METADATA = SCRIPT_DIR / 'ncbi_dataset.tsv'

df = pd.read_csv(NCBI_METADATA, sep='\t')
df['Organism Name'] = df['Organism Name'].str.replace(r'[^a-zA-Z0-9]', '_', regex=True)

if df['Organism Name'].is_unique:
    print("Organism Name is unique")
else:
    print("Organism Name is not unique")
    print(df['Organism Name'][df['Organism Name'].duplicated()].unique())

df[['Assembly Accession', 'Organism Name']].to_csv(SCRIPT_DIR / 'ncbi_dataset_only_organism.tsv', index=False)
