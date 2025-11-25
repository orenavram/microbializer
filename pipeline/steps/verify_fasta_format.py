import argparse
from pathlib import Path
import sys
import re
import collections
from Bio import SeqIO
import Bio.SeqUtils

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries import consts


ILLEGAL_CHARS_IN_RECORD_IDS = ':;,\'\"'


def record_has_illegal_chars(record_header):
    # Check if the first word in the string (which is the record ID) contains illegal characters
    record_id = re.match(r">(\S+)", record_header)
    if record_id is None:
        return True

    record_id = record_id.group(1)
    return any(not char.isascii() or char in ILLEGAL_CHARS_IN_RECORD_IDS for char in record_id)


def verify_fasta_format(logger, file_path, inputs_fasta_type):
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(
        Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)

    if file_path.suffix in ['.zip', '.gz', '.rar']:
        return f'{file_path.name} is a binary file (rather than textual). Please upload your genomes as FASTA text ' \
               f'files (such as fas or fna).'
    strain_name = file_path.stem
    with open(file_path) as f:
        line_number = 0
        try:
            line = f.readline()
            line_number += 1
            if not line:
                return f'Illegal FASTA format. First line in "{file_path.name}" is empty.'
            if not line.startswith('>'):
                return f'Illegal FASTA format. First line in "{file_path.name}" starts with "{line[0]}" instead of ">".'
            if record_has_illegal_chars(line):
                return f'Illegal format. First line in "{file_path.name}" contains an illegal character in its ' \
                       f'first word (one of: {" ".join(ILLEGAL_CHARS_IN_RECORD_IDS)}).'

            curated_content = f'>{strain_name}:{line[1:]}'
            previous_line_was_header = True

            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    if previous_line_was_header:
                        return f'Illegal FASTA format. "{file_path.name}" contains an empty record. ' \
                               f'Both lines {line_number - 1} and {line_number} start with ">".'
                    elif record_has_illegal_chars(line):
                        return f'Illegal format. Line {line_number} in "{file_path.name}" contains an illegal ' \
                               f'character in its first word (one of: {" ".join(ILLEGAL_CHARS_IN_RECORD_IDS)}).'
                    else:
                        curated_content += f'>{strain_name}:{line[1:]}\n'
                        previous_line_was_header = True
                        continue

                else:  # not a header
                    previous_line_was_header = False
                    for c in line:
                        if c not in legal_chars:
                            return f'Illegal FASTA format. Line {line_number} in "{file_path.name}" contains ' \
                                   f'illegal DNA character "{c}".'
                    curated_content += f'{line}\n'

        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'Illegal FASTA format. Line {line_number} in "{file_path.name}" contains one (or more) non ' \
                   f'ascii character(s).'

    # Now that we verified the fasta format, we parse it again with Bio.SeqIO.
    record_ids = []
    total_genome_length = 0
    max_record_length = 0
    max_record_length_name = ''
    for record in SeqIO.parse(file_path, 'fasta'):
        record_ids.append(record.id)
        total_genome_length += len(record.seq)
        if len(record.seq) > max_record_length:
            max_record_length = len(record.seq)
            max_record_length_name = record.id

    duplicate_ids = [item for item, count in collections.Counter(record_ids).items() if count > 1]
    if duplicate_ids:
        return f'Illegal FASTA format. "{file_path.name}" contains duplicated record ids: {",".join(duplicate_ids)}.'

    if total_genome_length < consts.MIN_GENOME_LENGTH and inputs_fasta_type == 'genomes':
        return (
            f'Each FASTA file should contain the genome of a bacterium, hence it must contain at least {consts.MIN_GENOME_LENGTH} '
            f'nucleotides. It is a requirement for Prodigal to run successfully (the first step of the pipeline '
            f'that predicts ORFs from the genomes). {file_path.name} contains less than {consts.MIN_GENOME_LENGTH} nucleotides.')

    if max_record_length > consts.MAX_ORF_LENGTH and inputs_fasta_type == 'orfs':
        return (f'The ORF {max_record_length_name} in the input FASTA file {file_path.name} is longer than '
                f'{consts.MAX_ORF_LENGTH} nucleotides. This is biologically invalid since the longest known '
                f'bacterial or archeal ORF is {consts.MAX_ORF_LENGTH} nucleotides long.')

    # override the old file with the curated content
    with open(file_path, 'w') as f:
        f.write(curated_content)


def main(logger, job_input_path, errors_dir, inputs_fasta_type):
    with open(job_input_path, 'r') as f:
        for line in f:
            file_path = Path(line.strip())
            verification_error = verify_fasta_format(logger, file_path, inputs_fasta_type)
            if verification_error:
                with open(errors_dir / f'{file_path.stem}_error.txt', 'w') as error_f:
                    error_f.write(verification_error)
                break


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('errors_dir', type=Path, help='path to the errors in fasta format dir')
    parser.add_argument('inputs_fasta_type', type=str)
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, main, args.job_input_path, args.errors_dir, args.inputs_fasta_type)
