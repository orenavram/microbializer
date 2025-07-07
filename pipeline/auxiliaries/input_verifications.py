import Bio.SeqUtils
import shutil
import tarfile
import collections
import sys
import subprocess
from pathlib import Path
import re
from Bio import SeqIO

from . import consts
from .general_utils import remove_path, fail, write_done_file
from .configuration import Config

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.flask.flask_interface_consts import WEBSERVER_NAME

ILLEGAL_CHARS_IN_FILE_NAMES = '\\/:*?\"\'<>` |(),;&'
ILLEGAL_CHARS_IN_RECORD_IDS = ':;,\'\"'


def record_has_illegal_chars(record_header):
    # Check if the first word in the string (which is the record ID) contains illegal characters
    record_id = re.match(r">(\S+)", record_header)
    if record_id is None:
        return True

    record_id = record_id.group(1)
    return any(not char.isascii() or char in ILLEGAL_CHARS_IN_RECORD_IDS for char in record_id)


def prepare_and_verify_input_data(logger, config: Config):
    done_file_path = config.done_files_dir / '00_prepare_and_verify_inputs.txt'
    if not done_file_path.exists():
        logger.info(f'Examining {config.raw_data_path}..')

        # extract zip and detect data folder
        primary_data_path = unpack_data(logger, config.raw_data_path, config.run_dir, config.error_file_path)

        remove_system_files(logger, primary_data_path)

        # copies input contigs_dir because we edit the files and want to keep the input directory as is
        shutil.copytree(primary_data_path, config.data_path, dirs_exist_ok=True)
        logger.info(f'Copied {primary_data_path} to {config.data_path}')

        # have to be AFTER system files removal (in the weird case a file name starts with a space)
        genome_names = set()
        for file_path in config.data_path.iterdir():
            new_file_path = fix_file_path(logger, config, file_path)
            genome_name = new_file_path.stem
            if genome_name in genome_names:
                error_msg = f'Two (or more) of the uploaded genomes has the same name ' \
                            f'e.g., {genome_name}. Please make sure each file name is unique.'
                fail(logger, error_msg, config.error_file_path)
            genome_names.add(genome_name)

        data_path_files = list(config.data_path.iterdir())
        number_of_genomes = len(data_path_files)
        logger.info(f'Number of genomes to analyze is {number_of_genomes}')
        logger.info(f'data_path contains the following {number_of_genomes} files: {data_path_files}')

        # check MINimal number of genomes
        if number_of_genomes < consts.MIN_NUMBER_OF_GENOMES_TO_ANALYZE and not config.bypass_number_of_genomes_limit:
            error_msg = f'The dataset contains too few genomes ({WEBSERVER_NAME} does comparative analysis and ' \
                        f'thus needs at least 2 genomes).'
            fail(logger, error_msg, config.error_file_path)

        # check MAXimal number of genomes
        if number_of_genomes > consts.MAX_NUMBER_OF_GENOMES_TO_ANALYZE and not config.bypass_number_of_genomes_limit:
            error_msg = f'The dataset contains too many genomes. {WEBSERVER_NAME} allows analyzing up to ' \
                        f'{consts.MAX_NUMBER_OF_GENOMES_TO_ANALYZE} genomes due to the high resource consumption. However, ' \
                        f'upon request, we might be able to analyze larger datasets. Please contact us ' \
                        f'directly and we will try to do that for you.'
            fail(logger, error_msg, config.error_file_path)

        # must be only after the spaces removal from the species names!!
        verification_error = verify_fasta_format(logger, config.data_path, config.inputs_fasta_type)
        if verification_error:
            remove_path(logger, config.data_path)
            fail(logger, verification_error, config.error_file_path)

        genomes_names = [genome_path.stem for genome_path in sorted(config.data_path.iterdir())]

        with open(config.genomes_names_path, 'w') as genomes_name_fp:
            genomes_name_fp.write('\n'.join(genomes_names))

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def unpack_data(logger, raw_data_path, run_dir, error_file_path):
    if not raw_data_path.is_dir():
        unzipped_data_path = run_dir / 'data'
        try:
            if tarfile.is_tarfile(raw_data_path):
                logger.info(f'UnTARing {raw_data_path}...')
                with tarfile.open(raw_data_path, 'r:gz') as f:
                    f.extractall(path=unzipped_data_path)  # unzip tar folder to parent dir
                logger.info(f'Succeeded unTARing {raw_data_path} to {unzipped_data_path}')
            elif raw_data_path.suffix == '.gz':  # gunzip gz file
                cmd = f'gunzip -f "{raw_data_path}"'
                logger.info(f'Calling: {cmd}')
                subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                unzipped_data_path = raw_data_path.with_suffix("")  # trim the ".gz"
                logger.info(f'Succeeded gunzipping {raw_data_path} to {unzipped_data_path}')
            else:
                logger.info(f'UnZIPing {raw_data_path} using shutil.unpack_archive...')
                shutil.unpack_archive(raw_data_path, extract_dir=unzipped_data_path)  # unzip tar folder to parent dir
                logger.info(f'Succeeded unZIPing {raw_data_path} to {unzipped_data_path}')
        except Exception as e:
            logger.info(e)
            remove_path(logger, raw_data_path)
            fail(logger,
                 f'{WEBSERVER_NAME} failed to decompress your data. Please make sure all your FASTA files names '
                 f'do contain only dashes, dots, and alphanumeric characters (a-z, A-Z, 0-9). Other characters such as '
                 f'parenthesis, pipes, slashes, are not allowed. Please also make sure your archived file format is legal (either a '
                 f'.zip file or a .tar.gz file in which each file is in FASTA format containing genomic sequence of a different species).',
                 error_file_path)

        if not unzipped_data_path.exists():
            fail(logger, f'Failed to unzip {raw_data_path.name} (maybe it is empty?)', error_file_path)

        if not unzipped_data_path.is_dir():
            fail(logger, f'{raw_data_path.name} is not a zipped folder', error_file_path)

        remove_system_files(logger, unzipped_data_path)

        content = list(unzipped_data_path.iterdir())
        if len(content) == 1 and content[0].is_dir():
            data_path = content[0]
            logger.info(f'{unzipped_data_path} contains only 1 directory: {data_path.name}. Using it as data_path.')
        else:
            data_path = unzipped_data_path
    else:
        data_path = raw_data_path

    logger.info(f'Updated data_path is: {data_path}')
    remove_system_files(logger, data_path)

    if not list(data_path.iterdir()):
        fail(logger, f'No input files were found in {raw_data_path.name}', error_file_path)

    for file_path in data_path.iterdir():
        if file_path.is_dir():
            if file_path.name == '__MACOSX':
                # happens too many times to mac users so I decided to assist in this case
                logger.info('Removing __MACOSX dir...')
                shutil.rmtree(file_path, ignore_errors=True)
            else:
                fail(logger, f'Please make sure to upload one archived folder containing (only) FASTA files '
                             f'("{file_path.name}" is a folder).', error_file_path)
        else:
            # make sure each fasta has writing permissions for downstream editing
            file_path.chmod(0o644)  # -rw-r--r--

    return data_path


def fix_file_path(logger, config, file_path):
    file_name = file_path.name
    new_file_name = []
    for char in file_name:
        # Replace non-ASCII or illegal characters with _
        if not char.isascii() or char.isspace() or char in ILLEGAL_CHARS_IN_FILE_NAMES:
            new_file_name.append('_')
            logger.info(f'Illegal character ({char}) were detected in {file_name}! Replacing with "_".')
        else:
            new_file_name.append(char)

    new_file_name = ''.join(new_file_name)

    # Due to weird mmseqs behavior, the ORFs headers (and therefore the filenames) cannot start with 'consensus'.
    if new_file_name.startswith('consensus'):
        new_file_name = new_file_name.replace('consensus', 'consenzus')

    if file_name != new_file_name:
        # illegal character in file name were found
        new_file_path = rename_file(logger, file_path, new_file_name, config.error_file_path)
        if config.outgroup == file_path.stem:
            new_outgroup = Path(new_file_name).stem
            logger.info(f'Following the change of input genome name {file_name} to {new_file_name}, '
                        f'outgroup argument was changed from {config.outgroup} to {new_outgroup}')
            config.outgroup = new_outgroup
    else:
        new_file_path = file_path

    return new_file_path


def remove_system_files(logger, dir_path):
    for system_file_path in dir_path.glob("[._]*"):
        logger.warning(f'Removing system file: {system_file_path}')
        remove_path(logger, system_file_path)


def rename_file(logger, file_path, new_file_name, error_file_path):
    # a name replacement was applied. move the file to its new name.
    try:
        new_file_path = file_path.parent / new_file_name
        logger.info(f'Renaming file path from: {file_path.name} to {new_file_name}')
        file_path.rename(new_file_path)
        return new_file_path
    except:
        error_msg = f'One (or more) file name(s) contain illegal character such as parenthesis, pipes, or slashes.\n' \
                    f'In order to avoid downstream parsing errors, {WEBSERVER_NAME} ' \
                    f'automatically replaces these spaces with dashes. For some reason, the replacement for ' \
                    f'{file_path.name} failed. Please make sure all your input files names contain only dashes, dots, ' \
                    f'and alphanumeric characters (a-z, A-Z, 0-9) and re-submit your job.'
        fail(logger, error_msg, error_file_path)


def verify_fasta_format(logger, data_path, inputs_fasta_type):
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(
        Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    for file_path in data_path.iterdir():
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
            return (f'Each FASTA file should contain the genome of a bacterium, hence it must contain at least {consts.MIN_GENOME_LENGTH} '
                    f'nucleotides. It is a requirement for Prodigal to run successfully (the first step of the pipeline '
                    f'that predicts ORFs from the genomes). {file_path.name} contains less than {consts.MIN_GENOME_LENGTH} nucleotides.')

        if max_record_length > consts.MAX_ORF_LENGTH and inputs_fasta_type == 'orfs':
            return (f'The ORF {max_record_length_name} in the input FASTA file {file_path.name} is longer than '
                    f'{consts.MAX_ORF_LENGTH} nucleotides. This is biologically invalid since the longest known '
                    f'bacterial or archeal ORF is {consts.MAX_ORF_LENGTH} nucleotides long.')

        # override the old file with the curated content
        with open(file_path, 'w') as f:
            f.write(curated_content)
