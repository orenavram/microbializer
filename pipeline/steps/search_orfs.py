import subprocess
import sys
import argparse
from pathlib import Path
import mmap
import json

import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step


def find_genes(logger, genome_path, orfs_output_file_path):
    """
        input:path to fasta file with prokaryotic genome to be analyzed
        output: protein-coding gene prediction for input genome
    """
    cmd = f'prodigal -q -i "{genome_path}" -d "{orfs_output_file_path}" -o /dev/null'
    logger.info(f'Starting prodigal. Executed command is: {cmd}')
    subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

    if not orfs_output_file_path.exists() or orfs_output_file_path.stat().st_size == 0:
        raise Exception(f'Could not extract ORFs for {genome_path}')
    with open(orfs_output_file_path, 'rb', 0) as orf_f, mmap.mmap(orf_f.fileno(), 0, access=mmap.ACCESS_READ) as s:
        if s.find(b'>') == -1:
            raise Exception(f'{genome_path} does not contain any ORFs')

    logger.info(f'ORFs were extracted successfully for {genome_path} to {orfs_output_file_path}')


def mimic_prodigal_output(logger, input_orfs_path, output_orf_path):
    # edit headers of ORFs to match the structure of prodigal output
    fixed_content = ''
    with open(input_orfs_path, 'r') as orfs_file:
        for line in orfs_file:
            if line.startswith('>'):
                fixed_content += f'{line.strip()} # START # END # 1 # \n'
            else:
                fixed_content += line

    # override the old file with the fixed content
    with open(output_orf_path, 'w') as f:
        f.write(fixed_content)

    logger.info(f'Inputs are already orfs. Copied {input_orfs_path} to {output_orf_path}.')


def extract_orfs_statistics(logger, orf_path, orfs_statistics_dir):
    total_num_of_nucleotides = 0
    total_num_of_GC = 0
    orfs_count = 0
    for seq_record in SeqIO.parse(orf_path, 'fasta'):
        orfs_count += 1
        total_num_of_nucleotides += len(seq_record)
        total_num_of_GC += seq_record.count('G') + seq_record.count('C') + seq_record.count('g') + seq_record.count('c')

    orfs_statistics = {}
    orfs_statistics['orfs_count'] = orfs_count
    orfs_statistics['gc_content'] = (total_num_of_GC / total_num_of_nucleotides) * 100

    genome_name = orf_path.stem
    orfs_statistics_file_name = f'{genome_name}.json'
    orfs_statistics_file_path = orfs_statistics_dir / orfs_statistics_file_name
    with open(orfs_statistics_file_path, 'w') as fp:
        json.dump(orfs_statistics, fp)

    logger.info(f'ORFs statistics were extracted successfully for {genome_name} to {orfs_statistics_file_path}')


def fna_to_faa(logger, nucleotide_path, protein_path):
    with open(nucleotide_path, "r") as in_handle, open(protein_path, "w") as out_handle:
        # Iterate through each sequence record in the input file
        for record in SeqIO.parse(in_handle, "fasta"):
            # Translate the DNA sequence into a protein sequence
            translated_record = record.translate(id=True, name=True, description=True)

            # Write the translated record to the output file
            SeqIO.write(translated_record, out_handle, "fasta")

    logger.info(f'Translated fatsa file {nucleotide_path}. Output was written successfully to: {protein_path}')


def create_coordinates_file(logger, orfs_file_path, coordinates_file_path):
    record_ids = [record.id for record in SeqIO.parse(orfs_file_path, 'fasta')]
    coordinates = list(range(len(record_ids)))

    coordinates_df = pd.DataFrame({'record_id': record_ids, 'coordinate': coordinates})
    coordinates_df.to_csv(coordinates_file_path, index=False)
    logger.info(f'Coordinates file was created successfully for {orfs_file_path} in {coordinates_file_path}')


def find_genes_of_all_files(logger, job_input_path, orfs_sequences_dir, orfs_statistics_dir, orfs_translated_dir,
                            orfs_coordinates_dir, inputs_fasta_type):
    with open(job_input_path, 'r') as f:
        for line in f:
            genome_path = Path(line.strip())
            genome_name = genome_path.stem
            orfs_output_file_path = orfs_sequences_dir / f'{genome_name}.fna'

            if inputs_fasta_type == 'genomes':
                find_genes(logger, genome_path, orfs_output_file_path)
            else:  # inputs_fasta_type == 'orfs'
                mimic_prodigal_output(logger, genome_path, orfs_output_file_path)

            extract_orfs_statistics(logger, orfs_output_file_path, orfs_statistics_dir)
            fna_to_faa(logger, orfs_output_file_path, orfs_translated_dir / f'{genome_name}.faa')
            create_coordinates_file(logger, orfs_output_file_path, orfs_coordinates_dir / f'{genome_name}.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path,
                        help='path to a file that contains the genome names to search orfs for')
    parser.add_argument('orfs_sequences_dir', type=Path, help='path to orfs sequences dir')
    parser.add_argument('orfs_statistics_dir', type=Path, help='path to orfs statistics dir')
    parser.add_argument('orfs_translated_dir', type=Path, help='path to orfs translated dir')
    parser.add_argument('orfs_coordinates_dir', type=Path, help='path to orfs coordinates dir')
    parser.add_argument('inputs_fasta_type', help='')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, find_genes_of_all_files, args.job_input_path, args.orfs_sequences_dir, args.orfs_statistics_dir,
             args.orfs_translated_dir, args.orfs_coordinates_dir, args.inputs_fasta_type)
