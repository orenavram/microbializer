from pathlib import Path
import sys
import argparse
import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step


def get_genome_numeric_representation(logger, orthogroups_table_path, ORFs_dir_path, output_dir):
    orthogroups_df = pd.read_csv(orthogroups_table_path, index_col=0)

    genome_name_to_numeric_genome = {}
    for genome_name in orthogroups_df.columns:
        gene_id_to_og_number = {}
        for og, gene_ids in orthogroups_df[genome_name].items():
            if pd.isna(gene_ids):
                continue
            for gene_id in gene_ids.split(';'):
                gene_id = gene_id.strip()
                gene_id_to_og_number[gene_id] = og.lstrip('OG_')

        genome_orfs_path = ORFs_dir_path / f'{genome_name}.fna'
        orf_ids = [record.id for record in SeqIO.parse(genome_orfs_path, 'fasta')]
        numeric_genome = [gene_id_to_og_number.get(orf_id, '-1') for orf_id in orf_ids]
        genome_name_to_numeric_genome[genome_name] = ','.join(numeric_genome)

    output_path = output_dir / 'genome_numeric_representation.txt'
    with open(output_path, 'w') as f:
        for genome_name, numeric_genome in genome_name_to_numeric_genome.items():
            f.write(f'>{genome_name}\n{numeric_genome}\n')

    logger.info(f'Numeric genomes written to {output_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('orthogroups_table_path', type=Path, help='A path to an ortholog table')
    parser.add_argument('ORFs_dir_path', type=Path, help='A path to a ORF directory')
    parser.add_argument('output_dir', type=Path, help='A path to dir where the numeric genomes will be written')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, get_genome_numeric_representation, args.orthogroups_table_path, args.ORFs_dir_path, args.output_dir)
