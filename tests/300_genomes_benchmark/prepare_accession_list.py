from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

GENOMES_DIR = r'C:\Users\TalPNB22\Downloads\genomic datasets\salmonella_300\inputs'


def main():
    all_accessions = []
    for genome_path in Path(GENOMES_DIR).iterdir():
        if genome_path.is_file():
            accession_id = '_'.join(genome_path.stem.split('_')[:2])
            all_accessions.append(accession_id)

    with open(SCRIPT_DIR / 'accession_list.txt', 'w') as f:
        f.write('\n'.join(all_accessions))


if __name__ == '__main__':
    main()
