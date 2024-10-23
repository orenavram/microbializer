from Bio import SeqIO
import time
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

FASTA_FILE_PATH = '/groups/pupko/yairshimony/microbializer_runs/example_data/outputs7/02_ORFs/Ecoli_APEC_O1.02_ORFs'
FASTA_FILE_PATH = r"C:\Users\yairs\OneDrive\Documents\University\Master\Microbilaizer\example datasets\example input and output\M1CR0B1AL1Z3R_example_data_outputs\02_ORFs\Ecoli_APEC_O1.02_ORFs"


def load_header2sequences_dict(fasta_path, get_length=False, upper_sequence=False):
    header_to_sequence_dict = {}
    seq_length = 0

    with open(fasta_path) as f:
        header = f.readline().lstrip('>').rstrip()
        sequence = ''
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_length = len(sequence)
                if upper_sequence:
                    header_to_sequence_dict[header] = sequence.upper()
                else:
                    # leave untouched
                    header_to_sequence_dict[header] = sequence
                header = line.lstrip('>')
                sequence = ''
            else:
                sequence += line

        # don't forget last record!!
        if sequence != '':

            if upper_sequence:
                header_to_sequence_dict[header] = sequence.upper()
            else:
                # leave untouched
                header_to_sequence_dict[header] = sequence

    if get_length:
        return header_to_sequence_dict, seq_length
    else:
        return header_to_sequence_dict


def test_biopython():
    start = time.time()
    sequences_biopython = {}
    for seq_record in SeqIO.parse(FASTA_FILE_PATH, 'fasta'):
        sequences_biopython[seq_record.id] = seq_record.seq
    print(f'Biopython took {time.time() - start} seconds')

    start = time.time()
    sequences_manual = load_header2sequences_dict(FASTA_FILE_PATH)
    print(f'Manual parsing took {time.time() - start} seconds')


if __name__ == '__main__':
    test_biopython()


