import os
from sys import argv
import argparse
import logging
import sys
from Bio import SeqIO
import traceback
import shutil
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger
from auxiliaries.logic_auxiliaries import flatten
from auxiliaries import consts


def write_og_dna_sequences_file(og_name, og_members, gene_to_sequence_dict, og_dna_path):
    og_sequences = ''
    for gene_name in og_members:
        og_sequences += f'>{gene_name}\n{gene_to_sequence_dict[gene_name]}\n'

    with open(og_dna_path, 'w') as f:
        f.write(og_sequences)
    logger.info(f'Extracted dna sequences for {og_name} to {og_dna_path}')


def write_og_aa_sequences_file(og_name, og_members, protein_to_sequence_dict, og_aa_path):
    og_sequences = ''
    for protein_name in og_members:
        og_sequences += f'>{protein_name}\n{protein_to_sequence_dict[protein_name]}\n'

    with open(og_aa_path, 'w') as f:
        f.write(og_sequences)
    logger.info(f'Extracted protein sequences for {og_name} to {og_aa_path}')


def reconstruct_msa(logger, sequences_file_path, output_file_path):
    # cmd = f'mafft-linsi --maxiterate {maxiterate} --localpair {sequences_file_path} > {output_file_path}'

    # --auto Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i and FFT-NS-2, according to data size.
    # --amino/--nuc tells mafft that's an amino acid/nucleotide (respectively) msa. If you let it decide by itself, it
    # might wrong on small data sets as they might look like dna but they are NOT! e.g.,
    # [orenavr2@powerlogin-be2 test]$ cat /groups/pupko/orenavr2/igomeProfilingPipeline/experiments/test/analysis/motif_inference/17b_03/unaligned_sequences/17b_03_clusterRank_215_uniqueMembers_2_clusterSize_252.81.faa
    # >seq_235_lib_12_len_12_counts_126.40626975097965
    # CNTDVACAAPGN
    # >seq_1112_lib_C8C_len_10_counts_126.40626975097965
    # CTTACAPVNC
    records = list(SeqIO.parse(sequences_file_path, 'fasta'))
    if len(records) == 1:
        logger.info(f'Only one sequence in {sequences_file_path}. Copying it to {output_file_path}')
        shutil.copy(sequences_file_path, output_file_path)
    else:
        cmd = f'mafft --auto --amino --quiet {sequences_file_path} > {output_file_path}'
        logger.info(f'Starting MAFFT. Executed command is: {cmd}')
        subprocess.run(cmd, shell=True)
        logger.info(f'Finished MAFFT. Output was written to {output_file_path}')


def induce_sequence(logger, aligned_aa_seq, dna_seq):
    result = ''
    dna_i = 0
    for aa_i in range(len(aligned_aa_seq)):
        if aligned_aa_seq[aa_i] == '-':
            result += '-' * 3
        else:
            result += dna_seq[dna_i:dna_i + 3]
            dna_i += 3

    # TODO: remove this checkup
    if len(aligned_aa_seq) * 3 != len(result):
        logger.error('$' * 80)
        logger.error('len(aa_seq)*3 != len(result)')
        logger.error(f'{len(aligned_aa_seq) * 3} != {len(result)}')
    # # fill with trailing gaps so each induced dna sequence is of the same length
    # result += (len(aa_seq)*3-len(result))*'-'
    return result


def induce_msa(logger, og_members, gene_name_to_dna_sequence_dict, aa_msa_path, output_path):
    gene_name_to_aligned_aa_sequence = {record.id: record.seq for record in SeqIO.parse(aa_msa_path, 'fasta')}

    og_sequences = ''
    for gene_name in og_members:
        aligned_aa_sequence = gene_name_to_aligned_aa_sequence[gene_name]
        induced_dna_sequence = induce_sequence(logger, aligned_aa_sequence, gene_name_to_dna_sequence_dict[gene_name])
        og_sequences += f'>{gene_name}\n{induced_dna_sequence}\n'

    with open(output_path, 'w') as f:
        f.write(og_sequences)
    logger.info(f'Induced OG aa alignment to dna alignment. Output was written to {output_path}')


def extract_orfs(logger, all_orfs_path, all_proteins_path, orthogroups_file_path, start_og_number, end_og_number,
                 ogs_dna_output_dir, ogs_aa_output_dir, ogs_aa_aligned_output_dir, ogs_induced_dna_aligned_output_dir,
                 delimiter):
    gene_to_sequence_dict = {}
    for seq_record in SeqIO.parse(all_orfs_path, 'fasta'):
        gene_to_sequence_dict[seq_record.id] = seq_record.seq
    logger.info(f'Loaded all ({len(gene_to_sequence_dict)}) gene sequences into memory')

    protein_to_sequence_dict = {}
    for seq_record in SeqIO.parse(all_proteins_path, 'fasta'):
        protein_to_sequence_dict[seq_record.id] = seq_record.seq
    logger.info(f'Loaded all ({len(protein_to_sequence_dict)}) protein sequences into memory')

    with open(orthogroups_file_path) as f:
        header_line = f.readline()  # Skip the header line
        for i, line in enumerate(f):
            if i < start_og_number or i > end_og_number:
                continue
            line_tokens = line.strip().split(delimiter)
            og_name = line_tokens[0]
            og_members = flatten([strain_genes.split(';') for strain_genes in line_tokens[1:] if strain_genes])
            if not og_members:
                raise ValueError(f'Failed to extract any sequence for {og_name}.')

            logger.info(f'Extracting sequences for {og_name} ({len(og_members)} members)...')

            og_dna_path = os.path.join(ogs_dna_output_dir, f'{og_name}.fna')
            write_og_dna_sequences_file(og_name, og_members, gene_to_sequence_dict, og_dna_path)

            og_aa_path = os.path.join(ogs_aa_output_dir, f'{og_name}.faa')
            write_og_aa_sequences_file(og_name, og_members, protein_to_sequence_dict, og_aa_path)

            og_aligned_aa_path = os.path.join(ogs_aa_aligned_output_dir, f'{og_name}.faa')
            reconstruct_msa(logger, og_aa_path, og_aligned_aa_path)

            og_induced_dna_aligned_path = os.path.join(ogs_induced_dna_aligned_output_dir, f'{og_name}.fna')
            induce_msa(logger, og_members, gene_to_sequence_dict, og_aligned_aa_path, og_induced_dna_aligned_path)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('all_orfs_path', help='path to a file of all ORFs of all genomes')
    parser.add_argument('all_proteins_path', help='path to a file of all proteins of all genomes')
    parser.add_argument('orthogroups_file_path', help='path of the orthogroups file')
    parser.add_argument('start_og_number', type=int, help='OG number to start from')
    parser.add_argument('end_og_number', type=int, help='OG number to end at. Inclusive.')
    parser.add_argument('ogs_dna_output_dir', help='path to an output directory of ogs dna')
    parser.add_argument('ogs_aa_output_dir', help='path to an output directory of ogs aa')
    parser.add_argument('ogs_aa_aligned_output_dir', help='path to an output directory of ogs aligned aa')
    parser.add_argument('ogs_induced_dna_aligned_output_dir', help='path to an output directory of ogs induced aligned dna (codon alignment)')
    parser.add_argument('--delimiter', help='orthologs table delimiter', default=consts.CSV_DELIMITER)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        extract_orfs(logger, args.all_orfs_path, args.all_proteins_path, args.orthogroups_file_path, args.start_og_number,
                     args.end_og_number, args.ogs_dna_output_dir, args.ogs_aa_output_dir,
                     args.ogs_aa_aligned_output_dir, args.ogs_induced_dna_aligned_output_dir, args.delimiter)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
