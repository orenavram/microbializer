import shutil
import tarfile
import sys
import subprocess
from pathlib import Path
from collections import defaultdict

from . import consts
from .general_utils import remove_path, fail, write_done_file
from .configuration import Config
from .run_step_utils import wait_for_results, prepare_directories, submit_job, submit_batch

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.flask.flask_interface_consts import WEBSERVER_NAME


ILLEGAL_CHARS_IN_FILE_NAMES = '\\/:*?\"\'<>` |(),;&'


def prepare_and_verify_input_data(logger, times_logger, config: Config):
    done_file_path = config.done_files_dir / '00_prepare_and_verify_inputs.txt'
    if not done_file_path.exists():
        logger.info(f'Examining {config.raw_data_path}..')

        # extract zip and detect data folder
        primary_data_path = unpack_data(logger, config.raw_data_path, config.run_dir, config.error_file_path)

        remove_system_files(logger, primary_data_path)

        # copies input contigs_dir because we edit the files and want to keep the input directory as is
        copy_dir_cmd = f'rsync -a {primary_data_path}/ {config.data_path}/'
        logger.info(f'Running: {copy_dir_cmd}')
        subprocess.run(copy_dir_cmd, shell=True, check=True, capture_output=True, text=True)
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
        verification_error = verify_fasta_format(logger, times_logger, config)
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


def verify_fasta_format(logger, times_logger, config):
    step_number = '00_0'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_verify_fasta_format'
    script_path = consts.SRC_DIR / 'steps' / 'verify_fasta_format.py'
    errors_dir, verify_fasta_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir, step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Verify fasta formats...')

        job_index_to_fasta_files = defaultdict(list)

        for i, fasta_file in enumerate(config.data_path.iterdir()):
            job_index = i % config.max_parallel_jobs
            job_index_to_fasta_files[job_index].append(str(fasta_file))

        jobs_inputs_dir = verify_fasta_tmp_dir / consts.STEP_INPUTS_DIR_NAME
        jobs_inputs_dir.mkdir(parents=True, exist_ok=True)

        for job_index, job_fasta_files in job_index_to_fasta_files.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_fasta_files))

        script_params = [errors_dir, config.inputs_fasta_type]
        submit_batch(logger, config, script_path, script_params, verify_fasta_tmp_dir, 'verify_fasta')

        wait_for_results(logger, times_logger, step_name, verify_fasta_tmp_dir, config.error_file_path)

        for error_file in errors_dir.iterdir():
            if error_file.is_file() and error_file.suffix == '.txt':
                with open(error_file, 'r') as ef:
                    error_message = ef.read()
                return error_message

        write_done_file(logger, done_file_path)
