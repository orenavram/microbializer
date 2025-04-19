import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

FULL_METADATA_PATH = r'C:\repos\microbializer\tests\6000_benchmark\ncbi_dataset_only_organism.tsv'
ORFS_DIR = SCRIPT_DIR / 'orfs'


def main():
    df = pd.read_csv(FULL_METADATA_PATH)
    organism_list = [file_path.stem for file_path in ORFS_DIR.iterdir()]
    df_filter = df[df['Organism Name'].isin(organism_list)]
    df_filter.to_csv(SCRIPT_DIR / 'accession_list.csv', index=False)


if __name__ == '__main__':
    main()
