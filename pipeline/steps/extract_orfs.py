import argparse
from pathlib import Path
import sys
from Bio import SeqIO
import shutil
import subprocess
from collections import Counter
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.general_utils import flatten, none_or_path


def write_og_dna_sequences_file(logger, og_name, og_members, gene_to_sequence_dict, og_dna_path):
    og_sequences = ''
    for gene_name in og_members:
        og_sequences += f'>{gene_name}\n{gene_to_sequence_dict[gene_name]}\n'

    with open(og_dna_path, 'w') as f:
        f.write(og_sequences)
    logger.info(f'Extracted dna sequences for {og_name} to {og_dna_path}')


def write_og_aa_sequences_file(logger, og_name, og_members, protein_to_sequence_dict, og_aa_path):
    og_sequences = ''
    for protein_name in og_members:
        og_sequences += f'>{protein_name}\n{protein_to_sequence_dict[protein_name]}\n'

    with open(og_aa_path, 'w') as f:
        f.write(og_sequences)
    logger.info(f'Extracted protein sequences for {og_name} to {og_aa_path}')


def reconstruct_msa(logger, sequences_file_path, output_file_path):
    # --auto Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i and FFT-NS-2, according to data size.
    # --amino/--nuc tells mafft that's an amino acid/nucleotide (respectively) msa. If you let it decide by itself, it
    # might wrong on small data sets as they might look like dna but they are NOT!
    records = list(SeqIO.parse(sequences_file_path, 'fasta'))
    max_record_length = max(len(record.seq) for record in records)
    if len(records) == 1:
        logger.info(f'Only one sequence in {sequences_file_path}. Copying it to {output_file_path}')
        shutil.copy(sequences_file_path, output_file_path)
    else:
        if max_record_length < 5000:
            mode = '--auto'
        else:  # For large sequences, the --auto mode in MAFFT demands a lot of memory and it makes the job crash. So I use --parttree which uses less memory.
            mode = '--parttree'
        cmd = f'mafft {mode} --amino --quiet {sequences_file_path} > {output_file_path}'
        logger.info(f'Starting MAFFT. Executed command is: {cmd}')
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        logger.info(f'Finished MAFFT. Output was written to {output_file_path}')


def calc_consensus(logger, og_aa_path, consensus_faa_path):
    # Read sequences into a list
    sequences = [str(record.seq) for record in SeqIO.parse(og_aa_path, "fasta")]
    number_of_sequences = len(sequences)
    alignment_length = len(sequences[0])  # All sequences should be the same length

    consensus = []
    for i in range(alignment_length):
        column = [seq[i] for seq in sequences]  # Extract the column (all bases at position i)

        counts = Counter(column)
        consensus_base, consensus_count = counts.most_common(1)[0]

        # Handle gaps (`-`)
        if consensus_base == '-':
            if consensus_count / number_of_sequences > 0.5:
                continue  # Skip column
            else:
                consensus_base = counts.most_common(2)[1][0]  # Find the most common base that is not a gap

        consensus.append(consensus_base)

    with open(consensus_faa_path, "w") as f:
        f.write(f">{og_aa_path.stem}_consensus\n{''.join(consensus)}\n")

    # og_tmp_dir = os.path.join(output_dir, f'{og_name}_tmp')
    # os.makedirs(og_tmp_dir, exist_ok=True)
    #
    # msa_stockholm_path = os.path.join(og_tmp_dir, f"{og_name}.sto")
    # AlignIO.convert(og_path, "fasta", msa_stockholm_path, "stockholm")
    #
    # msa_db_path = os.path.join(og_tmp_dir, f"{og_name}_msaDb")
    # profile_db_path = os.path.join(og_tmp_dir, f"{og_name}_profileDB")
    # consensus_db_path = os.path.join(og_tmp_dir, f"{og_name}_consensusDb")
    # consensus_faa_raw_path = os.path.join(og_tmp_dir, f"{og_name}_consensus.faa")
    #
    # cmds = [f"mmseqs convertmsa {msa_stockholm_path} {msa_db_path} -v 1",
    #         f"mmseqs msa2profile {msa_db_path} {profile_db_path} --match-mode 1 -v 1",
    #         f"mmseqs profile2consensus {profile_db_path} {consensus_db_path} -v 1",
    #         f"mmseqs result2flat {consensus_db_path} {consensus_db_path} {consensus_db_path} {consensus_faa_raw_path} -v 1"]
    #
    # for cmd in cmds:
    #     logger.info(f'Calling: {cmd}')
    #     subprocess.run(cmd, shell=True, check=True)
    #
    # record = SeqIO.parse(consensus_faa_raw_path, 'fasta').__next__()
    # record.id = f"{og_name}_consensus"
    #
    # SeqIO.write(record, consensus_faa_path, 'fasta')
    # shutil.rmtree(og_tmp_dir, ignore_errors=True)

    logger.info(f'Consensus calculation finished. Output written to {consensus_faa_path}')


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


def extract_orfs(logger, all_orfs_path, all_proteins_path, orthogroups_file_path, job_input_path,
                 ogs_dna_output_dir, ogs_aa_output_dir, ogs_aa_aligned_output_dir,
                 ogs_induced_dna_aligned_output_dir, ogs_aa_consensus):
    with open(job_input_path, 'r') as f:
        ogs_numbers = [line.strip() for line in f]

    if all_orfs_path:
        gene_to_sequence_dict = {}
        for seq_record in SeqIO.parse(all_orfs_path, 'fasta'):
            gene_to_sequence_dict[seq_record.id] = seq_record.seq
        logger.info(f'Loaded all ({len(gene_to_sequence_dict)}) gene sequences into memory')

    if all_proteins_path:
        protein_to_sequence_dict = {}
        for seq_record in SeqIO.parse(all_proteins_path, 'fasta'):
            protein_to_sequence_dict[seq_record.id] = seq_record.seq
        logger.info(f'Loaded all ({len(protein_to_sequence_dict)}) protein sequences into memory')

    orthogroups_df = pd.read_csv(orthogroups_file_path, dtype=str)
    orthogroups_df = orthogroups_df[orthogroups_df['OG_name'].isin(ogs_numbers)]

    for i, row in orthogroups_df.iterrows():
        og_name = row['OG_name']
        og_members = flatten([strain_genes.split(';') for strain_genes in row[1:].dropna()])
        if not og_members:
            raise ValueError(f'Failed to extract any sequence for {og_name}.')

        logger.info(f'Extracting sequences for {og_name} ({len(og_members)} members)...')

        if ogs_dna_output_dir:
            og_dna_path = ogs_dna_output_dir / f'{og_name}.fna'
            if not og_dna_path.exists():
                write_og_dna_sequences_file(logger, og_name, og_members, gene_to_sequence_dict, og_dna_path)

        if ogs_aa_output_dir:
            og_aa_path = ogs_aa_output_dir / f'{og_name}.faa'
            if not og_aa_path.exists():
                write_og_aa_sequences_file(logger, og_name, og_members, protein_to_sequence_dict, og_aa_path)

        if ogs_aa_aligned_output_dir:
            og_aligned_aa_path = ogs_aa_aligned_output_dir / f'{og_name}.faa'
            if not og_aligned_aa_path.exists():
                reconstruct_msa(logger, og_aa_path, og_aligned_aa_path)

        if ogs_induced_dna_aligned_output_dir:
            og_induced_dna_aligned_path = ogs_induced_dna_aligned_output_dir / f'{og_name}.fna'
            if not og_induced_dna_aligned_path.exists():
                induce_msa(logger, og_members, gene_to_sequence_dict, og_aligned_aa_path, og_induced_dna_aligned_path)

        if ogs_aa_consensus:
            og_aa_consensus_path = ogs_aa_consensus / f'{og_name}.faa'
            if not og_aa_consensus_path.exists():
                calc_consensus(logger, og_aligned_aa_path, og_aa_consensus_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('all_orfs_path', type=none_or_path, help='path to a file of all ORFs of all genomes')
    parser.add_argument('all_proteins_path', type=none_or_path, help='path to a file of all proteins of all genomes')
    parser.add_argument('orthogroups_file_path', type=Path, help='path of the orthogroups file')
    parser.add_argument('ogs_dna_output_dir', type=none_or_path, help='path to an output directory of ogs dna')
    parser.add_argument('ogs_aa_output_dir', type=none_or_path, help='path to an output directory of ogs aa')
    parser.add_argument('ogs_aa_aligned_output_dir', type=none_or_path,
                        help='path to an output directory of ogs aligned aa')
    parser.add_argument('ogs_induced_dna_aligned_output_dir', type=none_or_path,
                        help='path to an output directory of ogs induced aligned dna (codon alignment)')
    parser.add_argument('ogs_aa_consensus', type=none_or_path,
                        help='path to an output directory of ogs consensus aa sequence')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, extract_orfs, args.all_orfs_path, args.all_proteins_path, args.orthogroups_file_path,
             args.job_input_path, args.ogs_dna_output_dir, args.ogs_aa_output_dir, args.ogs_aa_aligned_output_dir,
             args.ogs_induced_dna_aligned_output_dir, args.ogs_aa_consensus)
