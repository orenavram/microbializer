import os
import shutil
import subprocess
from time import time, sleep
import re
import logging
from datetime import timedelta, datetime
import pandas as pd
import stat
import sys
import argparse
from pathlib import Path

from . import consts
from .email_sender import send_email
from .q_submitter_power import submit_cmds_from_file_to_q

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from flask import flask_interface_consts, SharedConsts


def validate_slurm_error_logs(logger, slurm_logs_dir, error_file_path):
    for slurm_log_file in slurm_logs_dir.glob('*.err'):
        if slurm_log_file.stat().st_size == 0:
            continue

        with open(slurm_log_file) as f:
            for line in f:
                pass
            last_line = line

        if last_line.startswith('slurmstepd: error'):
            fail(logger, f'file {slurm_log_file} shows a slurm error: {last_line}', error_file_path)


def wait_for_results(logger, times_logger, script_name, path, num_of_expected_results, error_file_path,
                     start=0, recursive_step=False):
    """waits until path contains num_of_expected_results .done files"""
    if not start:
        start = time()
    logger.info(f'Waiting for {script_name}... Continues when {num_of_expected_results} results will be in: {path}')

    if num_of_expected_results == 0:
        raise ValueError('Number of expected results is 0! Something went wrong in the previous analysis steps.')

    total_time = 0
    i = 0
    current_num_of_results = 0
    while num_of_expected_results > current_num_of_results:
        assert not error_file_path.exists()

        current_num_of_results = sum(1 for _ in path.glob(f'*{consts.JOB_DONE_FILE_SUFFIX}'))
        jobs_left = num_of_expected_results - current_num_of_results
        sleep(consts.CHECK_JOB_DONE_INTERVAL_SECONDS)
        total_time += consts.CHECK_JOB_DONE_INTERVAL_SECONDS
        i += 1
        if i % 5 == 0:  # print status every 5 cycles of $time_to_wait
            logger.info(
                f'\t{timedelta(seconds=total_time)} have passed since started waiting ('
                f'{num_of_expected_results} - {current_num_of_results} = {jobs_left} more files are still missing)')

    end = time()
    total_time_waited = timedelta(seconds=int(end - start))
    logger.info(f'Done waiting for: {script_name} (took {total_time_waited}).')

    validate_slurm_error_logs(logger, path, error_file_path)

    if not recursive_step:
        walltime_sum, cpus_used_per_job, log_files_without_times, log_files_without_cpus = get_jobs_cummulative_time(
            logger, path)
        times_logger.info(f'Step {script_name} took {total_time_waited}. '
                          f'There were {num_of_expected_results} jobs and '
                          f'cumulatively they took {walltime_sum} wallclock time (times {cpus_used_per_job} cpus used '
                          f'per job = {walltime_sum * cpus_used_per_job} wallclock time). ' +
                          (f'Times are not complete since files {log_files_without_times} do not have times records. '
                           if log_files_without_times else '') +
                          (f'Cpus are not complete since files {log_files_without_cpus} do not have cpus records. '
                           if log_files_without_cpus else ''))
    else:
        sub_steps_total_cummulative_time, log_files_without_times, sub_steps_cummulative_times \
            = get_recursive_step_cummulative_times(path)
        for step_name, step_time in sub_steps_cummulative_times.items():
            times_logger.info(f'Sub-step {step_name} took cumulatively {step_time} wallclock time')

        times_logger.info(f'Step {script_name} took {total_time_waited}. There were {num_of_expected_results} jobs and '
                          f'cumulatively they took {sub_steps_total_cummulative_time} wallclock time.' +
                          (f'Times are not complete since files {log_files_without_times} do not have times records'
                           if log_files_without_times else ''))

    assert not error_file_path.exists()


def get_job_time_from_log_file(log_file_content, pattern_for_runtime, pattern_for_cpus):
    runtime_string_match = pattern_for_runtime.search(log_file_content)
    if not runtime_string_match:
        return None, None

    runtime_string = runtime_string_match.group(1)

    if '-' in runtime_string:
        days_string, time_string = runtime_string.split('-')
    else:
        days_string = 0
        time_string = runtime_string
    hours, minutes, seconds = map(float, time_string.split(':'))
    runtime = timedelta(days=int(days_string), hours=hours, minutes=minutes, seconds=seconds)

    cpus_string_match = pattern_for_cpus.search(log_file_content)
    if cpus_string_match:
        cpus_string = cpus_string_match.group(1)
        cpus = int(cpus_string)
    else:
        cpus = None

    return runtime, cpus


def get_jobs_cummulative_time(logger, path):
    pattern_for_walltime = re.compile(f'{consts.JOB_WALL_TIME_KEY}(.+) TimeLimit')
    pattern_for_cpus = re.compile(f'{consts.JOB_CPUS_KEY}(.+) NumTasks')

    walltime_sum = timedelta()
    cpus_used_in_jobs = []
    log_files_without_times = []
    log_files_without_cpus = []
    for log_file_path in path.glob('*log.txt'):
        with open(log_file_path, 'r') as log_file:
            content = log_file.read()

            walltime, cpus_used_per_job = get_job_time_from_log_file(content, pattern_for_walltime, pattern_for_cpus)
            if walltime is None:
                log_files_without_times.append(log_file_path)
                continue
            walltime_sum += walltime
            if cpus_used_per_job is None:
                log_files_without_cpus.append(log_file_path)
                continue
            cpus_used_in_jobs.append(cpus_used_per_job)

    cpus_used_in_jobs = set(cpus_used_in_jobs)
    # Sometimes when we submit jobs with 1 cpu, 2 are allocated and it's not an error
    if cpus_used_in_jobs == {1, 2} or cpus_used_in_jobs == {2}:
        cpus_used_in_jobs = {1}
    if len(cpus_used_in_jobs) != 1:  # The way that jobs are submitted ensures that each job in a step uses the same number of cpus
        logger.error(
            f'Not all jobs used the same number of cpus in path {path}. cpus_used_in_jobs = {cpus_used_in_jobs}. '
            f'Setting cpus_used_in_jobs to 1 to avoid errors.')
        cpus_used_in_jobs = {1}

    return walltime_sum, next(iter(cpus_used_in_jobs)), log_files_without_times, log_files_without_cpus


def get_recursive_step_cummulative_times(path):
    pattern = re.compile(f'Step (.+) took (.+)\. .* per job = (.+) wallclock time')

    log_files_without_times = []
    sub_steps_cummulative_times = {}
    for log_file_path in path.glob('*times_log.txt'):
        with open(log_file_path, 'r') as log_file:
            content = log_file.read()

        regex_search_result = pattern.finditer(content)
        if regex_search_result is None:
            log_files_without_times.append(log_file_path)
            continue

        for match in regex_search_result:
            step_name, step_running_time, step_cummulative_time = match.group(1), match.group(2), match.group(3)
            if step_name in sub_steps_cummulative_times:
                sub_steps_cummulative_times[step_name] += pd.Timedelta(step_cummulative_time)
            else:
                sub_steps_cummulative_times[step_name] = pd.Timedelta(step_cummulative_time)

    sub_steps_total_cummulative_time = sum(sub_steps_cummulative_times.values(), pd.Timedelta(0))
    return sub_steps_total_cummulative_time, log_files_without_times, sub_steps_cummulative_times


def prepare_directories(logger, outputs_dir_prefix, tmp_dir_prefix, dir_name):
    outputs_dir = outputs_dir_prefix / dir_name
    logger.info(f'Creating {outputs_dir}')
    outputs_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = tmp_dir_prefix / dir_name
    logger.info(f'Creating {tmp_dir}')
    tmp_dir.mkdir(parents=True, exist_ok=True)

    return outputs_dir, tmp_dir


def fail(logger, error_msg, error_file_path):
    logger.error(error_msg)
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise ValueError(error_msg)


def submit_mini_batch(logger, config, script_path, mini_batch_parameters_list, logs_dir, job_name, num_of_cpus=1,
                      memory=None, time_in_hours=None, command_to_run_before_script=None,
                      alternative_error_file=None):
    """
    :param script_path:
    :param mini_batch_parameters_list: a list of lists. each sublist corresponds to a single command and contain its parameters
    :param logs_dir:
    :param job_name:
    :param queue_name:
    :param verbose:
    :param num_of_cpus:
    :return: an example command to debug on the shell
    """

    shell_cmds_as_str = ''

    if config.use_job_manager:
        # shell_cmds_as_str += f'source ~/.bashrc{new_line_delimiter}'
        conda_sh_path = consts.CONDA_INSTALLATION_DIR / 'etc' / 'profile.d' / 'conda.sh'
        shell_cmds_as_str += f'source {conda_sh_path}\n'
        shell_cmds_as_str += f'conda activate {consts.CONDA_ENVIRONMENT_DIR}\n'
        shell_cmds_as_str += f'export PATH=$CONDA_PREFIX/bin:$PATH\n'

    # PREPARING RELEVANT COMMANDS
    if command_to_run_before_script:
        shell_cmds_as_str += f'{command_to_run_before_script}\n'

    error_file_path = alternative_error_file if alternative_error_file else config.error_file_path
    for params in mini_batch_parameters_list:
        shell_cmds_as_str += ' '.join(
            ['python', str(script_path), *[str(param) for param in params],
             f'-v {config.verbose}', f'--logs_dir {logs_dir}', f'--error_file_path {error_file_path}',
             f'--job_name {job_name}']) + '\n'

    # GENERATE DONE FILE
    shell_cmds_as_str += f'touch {logs_dir / (job_name + consts.JOB_DONE_FILE_SUFFIX)}\n'

    # WRITING CMDS FILE
    cmds_path = logs_dir / f'{job_name}.sh'
    with open(cmds_path, 'w') as f:
        f.write(shell_cmds_as_str)

    if config.use_job_manager:
        # Add execution permissions to cmds_path
        current_permissions = cmds_path.stat().st_mode
        cmds_path.chmod(current_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        if memory is None:
            memory = config.job_default_memory
        submit_cmds_from_file_to_q(logger, job_name, cmds_path, logs_dir, config.queue_name, str(num_of_cpus),
                                   config.account_name, memory, time_in_hours, config.node_name)
    else:
        # fetch directly on shell
        for shell_cmd in shell_cmds_as_str.split('\n'):
            if shell_cmd:
                logger.info(f'Running command: {shell_cmd}')
                subprocess.run(shell_cmd, shell=True, check=True)


def submit_batch(logger, config, script_path, batch_parameters_list, logs_dir, job_name_suffix,
                 num_of_cmds_per_job=1, num_of_cpus=1, memory=None, time_in_hours=None):
    """
    :param script_path:
    :param batch_parameters_list: a list of lists. each sublist corresponds to a single command and contain its parameters
    :param logs_dir:
    :param job_name_suffix: a string that will be concatenated after the batch number as the job name
    :param queue_name:
    :param num_of_cmds_per_job:
    :param num_of_cpus:
    :return: number of batches submitted (in case waiting for the results) and an example command to debug on the shell
    """
    num_of_mini_batches = 0

    job_name_suffix = job_name_suffix.replace(' ', '_')  # job name cannot contain spaces!

    for i in range(0, len(batch_parameters_list), num_of_cmds_per_job):
        mini_batch_parameters_list = batch_parameters_list[i: i + num_of_cmds_per_job]
        mini_batch_job_name = f'{num_of_mini_batches}_{job_name_suffix}'
        submit_mini_batch(
            logger, config, script_path, mini_batch_parameters_list, logs_dir, mini_batch_job_name,
            num_of_cpus=num_of_cpus, memory=memory, time_in_hours=time_in_hours)

        num_of_mini_batches += 1
        sleep(0.1)

    return num_of_mini_batches


def send_email_in_pipeline_end(logger, process_id, email_address, job_name, state):
    email_addresses = [SharedConsts.OWNER_EMAIL]
    email_addresses.extend(flask_interface_consts.ADDITIONAL_OWNER_EMAILS)
    if email_address:
        email_addresses.append(email_address)
    else:
        logger.warning(f'process_id = {process_id} email_address is empty, state = {state}, job_name = {job_name}')

    # sends mail once the job finished or crashes
    if state == SharedConsts.State.Finished:
        send_email(logger, SharedConsts.EMAIL_CONSTS.create_title(state, job_name),
                   SharedConsts.EMAIL_CONSTS.CONTENT_PROCESS_FINISHED.format(process_id=process_id), email_addresses)
    elif state == SharedConsts.State.Crashed:
        send_email(logger, SharedConsts.EMAIL_CONSTS.create_title(state, job_name),
                   SharedConsts.EMAIL_CONSTS.CONTENT_PROCESS_CRASHED.format(process_id=process_id), email_addresses)


def add_results_to_final_dir(logger, source_dir_path, final_output_dir):
    dest = final_output_dir / consts.OUTPUTS_DIRECTORIES_MAP[source_dir_path.name]

    try:
        logger.info(f'Copying {source_dir_path} TO {dest}')
        shutil.copytree(source_dir_path, dest, dirs_exist_ok=True)
    except Exception:
        logger.exception(f'Failed to copy {source_dir_path} to {dest}')
        raise

    return dest


def remove_path(logger, path_to_remove):
    logger.info(f'Removing {path_to_remove} ...')
    try:
        shutil.rmtree(path_to_remove)  # maybe it's a folder
    except:
        pass
    try:
        os.remove(path_to_remove)
    except:
        pass


def get_job_logger(log_file_dir, job_name, verbose):
    job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, '')

    logger = logging.getLogger('main')
    file_handler = logging.FileHandler(log_file_dir / f'{job_name}_{job_id}_log.txt', mode='a')
    formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    return logger


def get_job_times_logger(log_file_dir, job_name, verbose):
    job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, '')

    logger = logging.getLogger('times')
    file_handler = logging.FileHandler(log_file_dir / f'{job_name}_{job_id}_times_log.txt', mode='a')
    formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    return logger


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ("yes", "true", "t", "1"):
        return True
    elif value.lower() in ("no", "false", "f", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def none_or_path(value):
    if value == 'None':
        return None
    return Path(value)


def submit_clean_folders_job(logger, config):
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

    logger.info('Submitting a job to clean old user results...')
    tmp_dir = consts.CLEAN_JOBS_LOGS_DIR / now.strftime(datetime_format)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    clean_old_jobs_error_file_path = tmp_dir / 'error.txt'
    script_path = consts.SRC_DIR / 'steps' / 'clean_old_jobs.py'
    params = []

    submit_mini_batch(logger, config, script_path, [params], tmp_dir, 'clean_old_jobs',
                      alternative_error_file=clean_old_jobs_error_file_path)


def add_default_step_args(args_parser):
    args_parser.add_argument('-v', '--verbose', help='Increase output verbosity', type=str_to_bool)
    args_parser.add_argument('--logs_dir', type=Path, help='path to tmp dir to write logs to')
    args_parser.add_argument('--error_file_path', type=Path, help='path to error file')
    args_parser.add_argument('--job_name', help='job name')


def write_done_file(logger, done_file_path):
    with open(done_file_path, 'w') as f:
        f.write('.')
    logger.info(f'{done_file_path} was generated.')
