import os
import Bio.SeqUtils
import logging
logger = logging.getLogger('main')

def verify_fasta_format(data_path):
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    for file_name in os.listdir(data_path):
        if file_name.endswith('zip') or file_name.endswith('gz') or file_name.endswith('rar'):
            return f'{file_name} is a binary file (rather than textual). Please upload your genomes as text files (such as fas or fna).'
        file_path = os.path.join(data_path, file_name)
        strain_name = os.path.splitext(file_name)[0]
        with open(file_path) as f:
            line_number = 0
            try:
                line = f.readline()  # TODO .replace(",", "_")
                line_number += 1
                if not line:
                    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. First line in "{file_name}" is empty.'
                if not line.startswith('>'):
                    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. First line in "{file_name}" starts with "{line[0]}" instead of ">".'
                previous_line_was_header = True
                putative_end_of_file = False
                curated_content = f'>{strain_name}_{line[1:]}'.replace("|", "_")
                for line in f:
                    line_number += 1
                    line = line.strip()
                    if not line:
                        # if not putative_end_of_file: # ignore trailing empty lines
                        #     putative_end_of_file = line_number
                        continue
                    if putative_end_of_file:  # non empty line after empty line
                        return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {putative_end_of_file} in "{file_name}" is empty.'
                    if line.startswith('>'):
                        if previous_line_was_header:
                            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. "{file_name}" contains an empty record. Both lines {line_number-1} and {line_number} start with ">".'
                        else:
                            previous_line_was_header = True
                            curated_content += f'>{strain_name}_{line[1:]}\n'.replace("|", "_")
                            continue
                    else:  # not a header
                        previous_line_was_header = False
                        for c in line:
                            if c not in legal_chars:
                                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in "{file_name}" contains illegal DNA character "{c}".'
                        curated_content += f'{line}\n'
            except UnicodeDecodeError as e:
                logger.info(e.args)
                line_number += 1  # the line that was failed to be read
                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in "{file_name}" contains one (or more) non <a href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s).'
        # override the old file with the curated content
        with open(file_path, 'w') as f:
            f.write(curated_content)