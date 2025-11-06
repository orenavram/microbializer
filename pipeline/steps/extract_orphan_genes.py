import argparse
from pathlib import Path
import sys
import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.general_utils import flatten


def extract_gene_names_from_fasta(proteins_file_path):
    with open(proteins_file_path) as proteins_file_fp:
        gene_names = [record.name.strip() for record in SeqIO.parse(proteins_file_fp, 'fasta')]
    return gene_names


def extract_orphan_proteins(logger, proteins_file_path, orthogroups_df, output_dir):
    strain_name = proteins_file_path.stem

    # Orphan ortogroups - orthogroups that contain genes from only one strain
    orphan_orthogroups = orthogroups_df[orthogroups_df.count(axis=1) == 1]
    orphan_orthogroups_of_strain = list(orphan_orthogroups[strain_name].dropna())
    genes_in_orphan_orthogroups_of_strain = flatten(
        [orthogroup.split(';') for orthogroup in orphan_orthogroups_of_strain])

    # Orphan genes (not in any orthogroup)
    orthogroups_strain_column = orthogroups_df[strain_name].dropna()
    genes_in_orthogroups = flatten([value.split(';') for value in orthogroups_strain_column])
    all_strain_genes = set(extract_gene_names_from_fasta(proteins_file_path))
    orphans = list(all_strain_genes.difference(genes_in_orthogroups))

    orphans_path = output_dir / f'{strain_name}_orphans.txt'
    with open(orphans_path, 'w') as orphans_path_fp:
        orphans_path_fp.write('\n'.join(orphan_orthogroups_of_strain + orphans))
    logger.info(f'Orphan genes of {strain_name} were written to {orphans_path}')

    orphans_count_path = output_dir / f'{strain_name}_orphans_stats.csv'
    orphans_stats = {
        'Orphan orthogroups count': len(orphan_orthogroups_of_strain),
        'Orphan single genes count': len(orphans),
        'Total orphans count': len(genes_in_orphan_orthogroups_of_strain) + len(orphans)
    }
    orphans_count_df = pd.DataFrame(orphans_stats, index=[strain_name])
    orphans_count_df.to_csv(orphans_count_path)
    logger.info(f'Orphan genes statistics of {strain_name} were written to {orphans_count_path}')


def extract_orphan_proteins_from_all_files(logger, job_input_path, orthogroups_file, output_dir):
    orthogroups_df = pd.read_csv(orthogroups_file)
    orthogroups_df.drop(columns=['OG_name'], inplace=True)

    with open(job_input_path, 'r') as f:
        for line in f:
            proteins_file_path = Path(line.strip())
            extract_orphan_proteins(logger, proteins_file_path, orthogroups_df, output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthogroups_file', type=Path, help='path to the the orthologs table')
    parser.add_argument('output_dir', type=Path, help='path to which the orphan proteins will be written')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, extract_orphan_proteins_from_all_files, args.job_input_path, args.orthogroups_file, args.output_dir)
