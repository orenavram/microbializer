import os
import Bio.SeqUtils
import shutil
import tarfile
import collections
import sys
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from . import consts
from .pipeline_auxiliaries import remove_path, fail
from .configuration import Config
from flask import flask_interface_consts


ILLEGAL_CHARS = '\\;:,^`~\'\"'


def has_illegal_chars(s):
    # Check if the first word in the string (which is the record ID) contains illegal characters
    record_id = s.split(' ')[0]
    return any(char in ILLEGAL_CHARS for char in record_id)


def prepare_and_verify_input_data(logger, config: Config):
    # extract zip and detect data folder
    primary_data_path = unpack_data(logger, config.raw_data_path, config.run_dir, config.error_file_path)

    for system_file_path in primary_data_path.glob("[._]*"):
        logger.warning(f'Removing system file: {system_file_path}')
        remove_path(logger, system_file_path)

    # copies input contigs_dir because we edit the files and want to keep the input directory as is
    shutil.copytree(primary_data_path, config.data_path, dirs_exist_ok=True)
    logger.info(f'Copied {primary_data_path} to {config.data_path}')

    # have to be AFTER system files removal (in the weird case a file name starts with a space)
    genome_names = set()
    for file_path in config.data_path.iterdir():
        new_file_name = fix_illegal_chars_in_file_name(logger, file_path.name)
        if file_path.name != new_file_name:
            # illegal character in file name were found
            move_file(logger, file_path, new_file_name, config.error_file_path)
            if config.outgroup == file_path.stem:
                new_outgroup = Path(new_file_name).stem
                logger.info(f'Following the change of input genome name {file_path.name} to {new_file_name}, '
                            f'outgroup argument was changed from {config.outgroup} to {new_outgroup}')
                config.outgroup = new_outgroup

        genome_name = file_path.stem
        if genome_name in genome_names:
            error_msg = f'Two (or more) of the uploaded geonmes contain the same name (prefix), ' \
                        f'e.g., {genome_name}. Please make sure each file name is unique.'
            fail(logger, error_msg, config.error_file_path)
        genome_names.add(genome_name)

    data_path_files = list(config.data_path.iterdir())
    number_of_genomes = len(data_path_files)
    logger.info(f'Number of genomes to analyze is {number_of_genomes}')
    logger.info(f'data_path contains the following {number_of_genomes} files: {data_path_files}')

    # check MINimal number of genomes
    if number_of_genomes < consts.MIN_NUMBER_OF_GENOMES_TO_ANALYZE and not config.bypass_number_of_genomes_limit:
        error_msg = f'The dataset contains too few genomes ({flask_interface_consts.WEBSERVER_NAME} does comparative analysis and ' \
                    f'thus needs at least 2 genomes).'
        fail(logger, error_msg, config.error_file_path)

    # check MAXimal number of genomes
    if number_of_genomes > consts.MAX_NUMBER_OF_GENOMES_TO_ANALYZE and not config.bypass_number_of_genomes_limit:
        error_msg = f'The dataset contains too many genomes. {flask_interface_consts.WEBSERVER_NAME} allows analyzing up to ' \
                    f'{consts.MAX_NUMBER_OF_GENOMES_TO_ANALYZE} genomes due to the high resource consumption. However, ' \
                    f'upon request, we might be able to analyze larger datasets. Please contact us ' \
                    f'directly and we will try to do that for you.'
        fail(logger, error_msg, config.error_file_path)

    # must be only after the spaces removal from the species names!!
    verification_error = verify_fasta_format(logger, config.data_path)
    if verification_error:
        remove_path(logger, config.data_path)
        fail(logger, verification_error, config.error_file_path)

    genomes_names = [genome_path.stem for genome_path in sorted(config.data_path.iterdir())]

    with open(config.genomes_names_path, 'w') as genomes_name_fp:
        genomes_name_fp.write('\n'.join(genomes_names))


def unpack_data(logger, raw_data_path, run_dir, error_file_path):
    if not raw_data_path.is_dir():
        unzipped_data_path = run_dir / 'data'
        try:
            if tarfile.is_tarfile(raw_data_path):
                logger.info('UnTARing')
                with tarfile.open(raw_data_path, 'r:gz') as f:
                    f.extractall(path=unzipped_data_path)  # unzip tar folder to parent dir
                logger.info('Succeeded!')
                # data_path = data_path.split('.tar')[0] # e.g., /groups/pupko/orenavr2/microbializer/example_data.tar.gz
                # logger.info(f'Updated data_path is:\n{data_path}')
            elif raw_data_path.suffix == '.gz':  # gunzip gz file
                cmd = f'gunzip -f "{raw_data_path}"'
                logger.info(f'Calling: {cmd}')
                subprocess.run(cmd, shell=True, check=True)
                unzipped_data_path = raw_data_path.with_suffix("")  # trim the ".gz"
            else:
                logger.info('UnZIPing')
                shutil.unpack_archive(raw_data_path, extract_dir=unzipped_data_path)  # unzip tar folder to parent dir
        except Exception as e:
            logger.info(e)
            remove_path(logger, raw_data_path)
            fail(logger, f'{flask_interface_consts.WEBSERVER_NAME} failed to decompress your data. Please make sure all your FASTA files names '
                 f'do contain only dashes, dots, and alphanumeric characters (a-z, A-Z, 0-9). Other characters such as '
                 f'parenthesis, pipes, slashes, are not allowed. Please also make sure your archived file format is legal (either a '
                 f'.zip file or a .tar.gz file in which each file is in FASTA format containing genomic sequence of a different species).',
                 error_file_path)
        logger.info('Succeeded!')

        if not unzipped_data_path.exists():
            fail(logger, f'Failed to unzip {raw_data_path.name} (maybe it is empty?)', error_file_path)

        if not unzipped_data_path.is_dir():
            fail(logger, 'Archived file content is not a folder', error_file_path)

        file = [x for x in os.listdir(unzipped_data_path) if not x.startswith(('_', '.'))][0]
        logger.info(f'first file in {unzipped_data_path} is:\n{file}')
        if (unzipped_data_path / file).is_dir():
            data_path = unzipped_data_path / file
            if not [x for x in os.listdir(data_path) if not x.startswith(('_', '.'))]:
                fail(logger, f'No input files were found in the archived folder.',
                     error_file_path)
        else:
            raw_data_path = unzipped_data_path

    logger.info(f'Updated data_path is: {raw_data_path}')
    for file_path in raw_data_path.iterdir():
        if file_path.suffix == '.gz':  # gunzip gz files in $data_path if any
            cmd = f'gunzip -f "{file_path}"'
            logger.info(f'Calling: {cmd}')
            subprocess.run(cmd, shell=True, check=True)

    for file_path in raw_data_path.iterdir():
        if not file_path.is_dir():
            # make sure each fasta has writing permissions for downstream editing
            os.chmod(file_path, 0o644)  # -rw-r--r--
        else:
            if file_path.name == '__MACOSX':
                # happens too many times to mac users so i decided to assist in this case
                logger.info('Removing __MACOSX file...')
                shutil.rmtree(file_path, ignore_errors=True)
            else:
                fail(logger, f'Please make sure to upload one archived folder containing (only) FASTA files '
                     f'("{file_path.name}" is a folder).', error_file_path)

    return raw_data_path


def fix_illegal_chars_in_file_name(logger, file_name, illegal_chars='\\|( )&:;,\xa0'):
    new_file_name = file_name
    for char in illegal_chars:
        if char in new_file_name:
            logger.info(f'File name with illegal character "{char}" was detected!\n')
            new_file_name = new_file_name.replace(char, '_')

    # Due to weird mmseqs behavior, the ORFs headers (and therefore the filenames) cannot start with 'consensus'.
    if new_file_name.startswith('consensus'):
        new_file_name = new_file_name.replace('consensus', 'consenzus')

    return new_file_name


def move_file(logger, file_path, new_file_name, error_file_path):
    # a name replacement was applied. move the file to its new name.
    try:
        new_file_path = file_path.parent / new_file_name
        logger.info(f'Renaming file path from:\n'
                    f'{file_path}\n'
                    f'to this:\n'
                    f'{new_file_path}')
        os.rename(file_path, new_file_path)
    except:
        error_msg = f'One (or more) file name(s) contain illegal character such as parenthesis, pipes, or slashes.\n' \
                    f'In order to avoid downstream parsing errors, {flask_interface_consts.WEBSERVER_NAME} ' \
                    f'automatically replaces these spaces with dashes. For some reason, the replacement for ' \
                    f'{file_path.name} failed. Please make sure all your input files names contain only dashes, dots, ' \
                    f'and alphanumeric characters (a-z, A-Z, 0-9) and re-submit your job.'
        fail(logger, error_msg, error_file_path)


def verify_fasta_format(logger, data_path):
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(
        Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    for file_path in data_path.iterdir():
        if file_path.suffix in ['.zip', '.gz', '.rar']:
            return f'{file_path.name} is a binary file (rather than textual). Please upload your genomes as text files ' \
                   f'(such as fas or fna).'
        strain_name = file_path.stem
        record_ids = []
        with open(file_path) as f:
            line_number = 0
            try:
                line = f.readline()  # TODO .replace(",", "_")
                line_number += 1
                if not line:
                    return f'Illegal FASTA format. First line in "{file_path.name}" is empty.'
                if not line.startswith('>'):
                    return f'Illegal FASTA format. First line in "{file_path.name}" starts with "{line[0]}" instead of ">".'
                if has_illegal_chars(line):
                    return f'Illegal format. First line in "{file_path.name}" contains an illegal character in its first word (one of: {ILLEGAL_CHARS}).'
                previous_line_was_header = True
                putative_end_of_file = False
                curated_content = f'>{strain_name}:{line[1:]}'
                record_ids.append(line[1:].split(' ')[0])
                for line in f:
                    line_number += 1
                    line = line.strip()
                    if not line:
                        # if not putative_end_of_file: # ignore trailing empty lines
                        #     putative_end_of_file = line_number
                        continue
                    if putative_end_of_file:  # non empty line after empty line
                        return f'Illegal FASTA format. Line {putative_end_of_file} in "{file_path.name}" is empty.'
                    if line.startswith('>'):
                        if previous_line_was_header:
                            return f'Illegal FASTA format. "{file_path.name}" contains an empty record. Both lines {line_number - 1} and {line_number} start with ">".'
                        elif has_illegal_chars(line):
                            return f'Illegal format. Line {line_number} in "{file_path.name}" contains an illegal character in its first word (one of: {ILLEGAL_CHARS}).'
                        else:
                            previous_line_was_header = True
                            curated_content += f'>{strain_name}:{line[1:]}\n'
                            record_ids.append(line[1:].split(' ')[0])
                            continue
                    else:  # not a header
                        previous_line_was_header = False
                        for c in line:
                            if c not in legal_chars:
                                return f'Illegal FASTA format. Line {line_number} in "{file_path.name}" contains illegal DNA character "{c}".'
                        curated_content += f'{line}\n'
            except UnicodeDecodeError as e:
                logger.info(e.args)
                line_number += 1  # the line that was failed to be read
                return f'Illegal FASTA format. Line {line_number} in "{file_path.name}" contains one (or more) non ascii character(s).'
        duplicate_ids = [item for item, count in collections.Counter(record_ids).items() if count > 1]
        if duplicate_ids:
            return f'Illegal FASTA format. "{file_path.name}" contains duplicated record ids: {",".join(duplicate_ids)}.'
        # override the old file with the curated content
        with open(file_path, 'w') as f:
            f.write(curated_content)
