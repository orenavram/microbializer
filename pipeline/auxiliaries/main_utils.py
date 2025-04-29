import shutil
from datetime import datetime, timedelta
import sys
from pathlib import Path
import pandas as pd
import math

from . import consts
from .run_step_utils import submit_mini_batch
from .email_sender import send_email

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.flask import flask_interface_consts, SharedConsts


def submit_clean_folders_job(logger, config):
    if not config.clean_intermediate_outputs:
        return

    logger.info('Cleaning up intermediate results...')

    folders_to_clean = [str(config.steps_results_dir)]

    clean_folders_tmp_dir = config.tmp_dir / 'clean_folders'
    clean_folders_tmp_dir.mkdir(parents=True, exist_ok=True)

    clean_folders_error_file_path = clean_folders_tmp_dir / 'error.txt'
    clean_folders_file_path = clean_folders_tmp_dir / 'folders.txt'
    with open(clean_folders_file_path, 'w') as fp:
        fp.write('\n'.join(folders_to_clean))

    params = [clean_folders_file_path]

    script_path = consts.SRC_DIR / 'steps' / 'clean_folders.py'
    submit_mini_batch(logger, config, script_path, [params], clean_folders_tmp_dir,
                      'clean_folders', alternative_error_file=clean_folders_error_file_path)


def submit_clean_old_user_results_job(logger, config):
    if not config.clean_old_job_directories:
        return

    logger.info('Checking if a clean old jobs job is needed...')
    consts.CLEAN_JOBS_LOGS_DIR.mkdir(parents=True, exist_ok=True)

    datetime_format = "%Y_%m_%d_%H_%M_%S"
    parsed_datetimes = []
    for subdir in consts.CLEAN_JOBS_LOGS_DIR.iterdir():
        if subdir.is_dir():
            try:
                parsed_datetimes.append(datetime.strptime(subdir.name, datetime_format))
            except ValueError:
                # Skip directories that don't match the datetime format
                pass

    now = datetime.now()
    if parsed_datetimes:
        most_recent = max(parsed_datetimes)
        if now - most_recent < timedelta(days=1):
            logger.info('No need to run a clean old user results job since the last run was less than a day ago.')
            return

    logger.info('Submitting a job to clean old user results since the last one was more than a day ago...')
    tmp_dir = consts.CLEAN_JOBS_LOGS_DIR / now.strftime(datetime_format)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    clean_old_jobs_error_file_path = tmp_dir / 'error.txt'
    script_path = consts.SRC_DIR / 'steps' / 'clean_old_jobs.py'
    params = []

    submit_mini_batch(logger, config, script_path, [params], tmp_dir, 'clean_old_jobs',
                      alternative_error_file=clean_old_jobs_error_file_path)


def send_email_in_pipeline_end(logger, config, state):
    if not config.send_email:
        logger.info('Not sending email since config.send_email is False')
        return

    email_addresses = [SharedConsts.OWNER_EMAIL]
    email_addresses.extend(flask_interface_consts.ADDITIONAL_OWNER_EMAILS)
    if config.email:
        email_addresses.append(config.email)
    else:
        logger.warning(f'process_id = {config.run_number} email_address is empty, state = {state}, job_name = {config.job_name}')

    logger.info(f'Sending email to {email_addresses} with state {state} and job_name {config.job_name}')

    # sends mail once the job finished or crashes
    if state.value == SharedConsts.State.Finished.value:
        send_email(logger, SharedConsts.EMAIL_CONSTS.create_title(state, config.job_name),
                   SharedConsts.EMAIL_CONSTS.CONTENT_PROCESS_FINISHED.format(process_id=config.run_number), email_addresses)
    elif state.value == SharedConsts.State.Crashed.value:
        send_email(logger, SharedConsts.EMAIL_CONSTS.create_title(state, config.job_name),
                   SharedConsts.EMAIL_CONSTS.CONTENT_PROCESS_CRASHED.format(process_id=config.run_number), email_addresses)


def add_results_to_final_dir(logger, source_dir_path, final_output_dir):
    dest = final_output_dir / consts.OUTPUTS_DIRECTORIES_MAP[source_dir_path.name]

    try:
        logger.info(f'Copying {source_dir_path} TO {dest}')
        shutil.copytree(source_dir_path, dest, dirs_exist_ok=True)
    except Exception:
        logger.exception(f'Failed to copy {source_dir_path} to {dest}')
        raise

    return dest


def initialize_progressbar(config):
    if config.only_calc_ogs:
        steps = consts.ONLY_CALC_OGS_TABLE_STEPS_NAMES_FOR_PROGRESS_BAR
    else:
        steps = consts.FULL_STEPS_NAMES_FOR_PROGRESS_BAR

    if not config.filter_out_plasmids:
        steps.remove('Filter out plasmids')

    if config.inputs_fasta_type == 'orfs':
        for i in range(len(steps)):
            if steps[i] == 'Predict and translate ORFs':
                steps[i] = 'Translate ORFs'
                break

    df = pd.DataFrame({'Step': steps, 'Finished': [False] * len(steps)})
    df.to_csv(config.progressbar_file_path, index=False)


def update_progressbar(progressbar_file_path, step_name_finished):
    df = pd.read_csv(progressbar_file_path)
    df.loc[df['Step'] == step_name_finished, 'Finished'] = True
    df.to_csv(progressbar_file_path, index=False)


def calc_genomes_batch_size(logger, config, num_of_genomes):
    if config.genomes_batch_size_calc_method == 'fixed_number':
        genomes_batch_size = config.genomes_batch_size
    elif config.genomes_batch_size_calc_method == 'sqrt':
        genomes_batch_size = math.ceil(math.sqrt(num_of_genomes))
    elif config.genomes_batch_size_calc_method == 'min_comparisons':
        genomes_batch_size = math.ceil((num_of_genomes * 2) ** (1/3))
    else:
        raise ValueError(f"Unknown genomes batch size calculation method: {config.genomes_batch_size_calc_method}")

    logger.info(f'Calculated genomes batch size: {genomes_batch_size}, according to method: '
                f'{config.genomes_batch_size_calc_method} and number of genomes: {num_of_genomes}')
    return genomes_batch_size


def define_intervals(number_of_genomes, genomes_batch_size):
    # Create the intervals
    intervals = []
    i = 0
    while i < number_of_genomes:
        interval_end = i + genomes_batch_size - 1
        intervals.append((i, interval_end))
        i = interval_end + 1

    # Adjust the last interval to ensure it ends exactly at end
    intervals[-1] = (intervals[-1][0], number_of_genomes - 1)

    return intervals


def zip_results(logger, config):
    if not config.step_to_complete and not config.only_calc_ogs and not config.do_not_copy_outputs_to_final_results_dir:
        logger.info('Zipping results folder...')
        shutil.make_archive(config.final_output_dir, 'zip', config.run_dir, config.final_output_dir_name)
        logger.info(f'Zipped results folder to {config.final_output_dir}.zip')
