import argparse
from pympler import asizeof

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('all_orfs_path', type=str, help='path to a file of all ORFs of all genomes')
    parser.add_argument('all_proteins_path', type=str, help='path to a file of all proteins of all genomes')
    args = parser.parse_args()

    gene_to_sequence_dict = {}
    for seq_record in SeqIO.parse(args.all_orfs_path, 'fasta'):
        gene_to_sequence_dict[seq_record.id] = seq_record.seq
    print(f'Loaded all ({len(gene_to_sequence_dict)}) gene sequences into memory')

    protein_to_sequence_dict = {}
    for seq_record in SeqIO.parse(args.all_proteins_path, 'fasta'):
        protein_to_sequence_dict[seq_record.id] = seq_record.seq
    print(f'Loaded all ({len(protein_to_sequence_dict)}) protein sequences into memory')

    genes_memory_mb = asizeof.asizeof(gene_to_sequence_dict) / (1024 ** 2)
    proteins_memory_mb = asizeof.asizeof(protein_to_sequence_dict) / (1024 ** 2)

    print(f'Gene sequences memory: {genes_memory_mb:.2f} MB')
    print(f'Protein sequences memory: {proteins_memory_mb:.2f} MB')


if __name__ == '__main__':
    main()
