import os
import shutil
import subprocess
from time import time, sleep
import re
import logging
from datetime import timedelta
import pandas as pd
import stat
import sys
import random
import string

from . import consts
from .email_sender import send_email
from .q_submitter_power import submit_cmds_from_file_to_q

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from flask import flask_interface_consts, SharedConsts


def execute(logger, process, process_is_string=False):
    process_str = process if process_is_string else ' '.join(str(token) for token in process)
    logger.info(f'Calling (process_is_string == {process_is_string}):\n{process_str}')
    subprocess.run(process, shell=process_is_string)


def validate_slurm_error_logs(logger, slurm_logs_dir, error_file_path):
    for file_name in os.listdir(slurm_logs_dir):
        if not file_name.endswith('.err'):
            continue

        file_path = os.path.join(slurm_logs_dir, file_name)
        if os.path.getsize(file_path) == 0:
            continue

        with open(file_path) as f:
            for line in f:
                pass
            last_line = line

        if last_line.startswith('slurmstepd: error'):
            fail(logger, f'file {file_path} shows a slurm error: {last_line}', error_file_path)


def wait_for_results(logger, times_logger, script_name, path, num_of_expected_results, error_file_path, suffix='done',
                     time_to_wait=10, start=0, error_message=None, recursive_step=False):
    """waits until path contains num_of_expected_results $suffix files"""
    if not start:
        start = time()
    logger.info(f'Waiting for {script_name}... Continues when {num_of_expected_results} results will be in: {path}')
    if num_of_expected_results == 0:
        if error_message:
            fail(logger, error_message, error_file_path)
        raise ValueError(
            f'\n{"#" * 100}\nNumber of expected results is {num_of_expected_results}! Something went wrong in the previous analysis steps...\n{"#" * 100}')
    total_time = 0
    i = 0
    current_num_of_results = 0
    while num_of_expected_results > current_num_of_results:
        assert not os.path.exists(error_file_path)
        try:
            current_num_of_results = sum(1 for x in os.listdir(path) if x.endswith(suffix))
        except:
            logger.info(f'Could not run the following command, probably due to some system error...')
            logger.info(f'current_num_of_results = sum(1 for x in os.listdir({path}) if x.endswith({suffix}))')
        jobs_left = num_of_expected_results - current_num_of_results
        sleep(time_to_wait)
        total_time += time_to_wait
        i += 1
        if i % 5 == 0:  # print status every 5 cycles of $time_to_wait
            logger.info(
                f'\t{timedelta(seconds=total_time)} have passed since started waiting ({num_of_expected_results} - {current_num_of_results} = {jobs_left} more files are still missing)')

    end = time()
    total_time_waited = timedelta(seconds=int(end - start))
    logger.info(f'Done waiting for: {script_name} (took {total_time_waited}).')

    validate_slurm_error_logs(logger, path, error_file_path)

    if not recursive_step:
        walltime_sum, cpus_used_per_job, log_files_without_times, log_files_without_cpus = get_jobs_cummulative_time(logger, path)
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

    assert not os.path.exists(error_file_path)


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
    log_files = [file_path for file_path in os.listdir(path) if file_path.endswith('log.txt')]
    pattern_for_walltime = re.compile(f'{consts.JOB_WALL_TIME_KEY}(.+) TimeLimit')
    pattern_for_cpus = re.compile(f'{consts.JOB_CPUS_KEY}(.+) NumTasks')

    walltime_sum = timedelta()
    cpus_used_in_jobs = []
    log_files_without_times = []
    log_files_without_cpus = []
    for log_file_name in log_files:
        log_file_full_path = os.path.join(path, log_file_name)
        with open(log_file_full_path, 'r') as log_file:
            content = log_file.read()

            walltime, cpus_used_per_job = get_job_time_from_log_file(content, pattern_for_walltime, pattern_for_cpus)
            if walltime is None:
                log_files_without_times.append(log_file_full_path)
                continue
            walltime_sum += walltime
            if cpus_used_per_job is None:
                log_files_without_cpus.append(log_file_full_path)
                continue
            cpus_used_in_jobs.append(cpus_used_per_job)

    cpus_used_in_jobs = set(cpus_used_in_jobs)
    if cpus_used_in_jobs == {1, 2} or cpus_used_in_jobs == {2}:  # Sometimes when we submit jobs with 1 cpu, 2 are allocated and it's not an error
        cpus_used_in_jobs = {1}
    if len(cpus_used_in_jobs) != 1:  # The way that jobs are submitted ensures that each job in a step uses the same number of cpus
        logger.error(f'Not all jobs used the same number of cpus in path {path}. cpus_used_in_jobs = {cpus_used_in_jobs}. '
                     f'Setting cpus_used_in_jobs to 1 to avoid errors.')
        cpus_used_in_jobs = {1}

    return walltime_sum, next(iter(cpus_used_in_jobs)), log_files_without_times, log_files_without_cpus


def get_recursive_step_cummulative_times(path):
    times_log_files = [file_path for file_path in os.listdir(path) if file_path.endswith('times_log.txt')]
    pattern = re.compile(f'Step (.+) took (.+)\. .* per job = (.+) wallclock time')

    log_files_without_times = []
    sub_steps_cummulative_times = {}
    for log_file_name in times_log_files:
        log_file_full_path = os.path.join(path, log_file_name)
        with open(log_file_full_path, 'r') as log_file:
            content = log_file.read()
        
        regex_search_result = pattern.finditer(content)
        if regex_search_result is None:
            log_files_without_times.append(log_file_full_path)
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
    outputs_dir = os.path.join(outputs_dir_prefix, dir_name)
    logger.info(f'Creating {outputs_dir}')
    os.makedirs(outputs_dir, exist_ok=True)

    tmp_dir = os.path.join(tmp_dir_prefix, dir_name)
    logger.info(f'Creating {tmp_dir}')
    os.makedirs(tmp_dir, exist_ok=True)

    return outputs_dir, tmp_dir


def fail(logger, error_msg, error_file_path):
    logger.error(error_msg)
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise ValueError(error_msg)


def submit_mini_batch(logger, script_path, mini_batch_parameters_list, logs_dir, error_file_path, queue_name,
                      account_name, job_name='', verbose=False, num_of_cpus=1, memory=consts.DEFAULT_MEMORY_PER_JOB_GB, time_in_hours=None,
                      done_file_is_needed=True, command_to_run_before_script=None, node_name=None):
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

    if consts.USE_JOB_MANAGER:
        # shell_cmds_as_str += f'source ~/.bashrc{new_line_delimiter}'
        conda_sh_path = os.path.join(consts.CONDA_INSTALLATION_DIR, 'etc/profile.d/conda.sh')
        shell_cmds_as_str += f'source {conda_sh_path}\n'
        shell_cmds_as_str += f'conda activate {consts.CONDA_ENVIRONMENT_DIR}\n'
        shell_cmds_as_str += f'export PATH=$CONDA_PREFIX/bin:$PATH\n'

    example_shell_cmd = ' '.join(['python', script_path, *[str(param) for param in mini_batch_parameters_list[0]]] + (
        ['-v'] if verbose else []))

    # PREPARING RELEVANT COMMANDS
    if command_to_run_before_script:
        shell_cmds_as_str += f'{command_to_run_before_script}\n'

    for params in mini_batch_parameters_list:
        shell_cmds_as_str += ' '.join(
            ['python', script_path, *[str(param) for param in params],
            '-v' if verbose else '', f'--logs_dir {logs_dir}', f'--error_file_path {error_file_path}']) + '\n'

    if not job_name:
        job_name = time()

    if done_file_is_needed:
        # GENERATE DONE FILE
        params = [os.path.join(logs_dir, job_name + '.done'), '', f'--logs_dir {logs_dir}']  # write an empty string (like "touch" command)
        shell_cmds_as_str += ' '.join(['python', os.path.join(consts.SRC_DIR, 'auxiliaries/file_writer.py'), *params]) + '\n'

    if consts.USE_JOB_MANAGER:
        # WRITING CMDS FILE
        cmds_path = os.path.join(logs_dir, f'{job_name}.sh')
        with open(cmds_path, 'w') as f:
            f.write(shell_cmds_as_str)

        # Add execution permissions to cmds_path
        current_permissions = os.stat(cmds_path).st_mode
        os.chmod(cmds_path, current_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        submit_cmds_from_file_to_q(logger, job_name, cmds_path, logs_dir, queue_name, str(num_of_cpus), account_name,
                                   memory, time_in_hours, node_name)
    else:
        # fetch directly on shell
        for shell_cmd in shell_cmds_as_str.split('\n'):
            execute(logger, shell_cmd, process_is_string=True)

    return example_shell_cmd


def submit_batch(logger, script_path, batch_parameters_list, logs_dir, error_file_path, job_name_suffix,
                 queue_name, account_name, num_of_cmds_per_job=1, num_of_cpus=1,
                 memory=consts.DEFAULT_MEMORY_PER_JOB_GB, time_in_hours=None, node_name=None):
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
    example_cmd_from_last_mini_batch = 'NO COMMANDS WERE FETCHED'

    job_name_suffix = job_name_suffix.replace(' ', '_')  # job name cannot contain spaces!

    for i in range(0, len(batch_parameters_list), num_of_cmds_per_job):
        mini_batch_parameters_list = batch_parameters_list[i: i + num_of_cmds_per_job]
        mini_batch_job_name = f'{num_of_mini_batches}_{job_name_suffix}'
        example_cmd_from_last_mini_batch = submit_mini_batch(logger, script_path, mini_batch_parameters_list, logs_dir,
                                                             error_file_path, queue_name, account_name, mini_batch_job_name,
                                                             verbose=False, num_of_cpus=num_of_cpus,
                                                             memory=memory, time_in_hours=time_in_hours, node_name=node_name)
        logger.info(f'Example command from batch {mini_batch_job_name}:\n{example_cmd_from_last_mini_batch}')
        num_of_mini_batches += 1
        sleep(0.1)

    return num_of_mini_batches, example_cmd_from_last_mini_batch


def wait_for_output_folder(logger, output_folder, max_waiting_time=300):
    i = 0
    while not os.path.exists(output_folder):
        logger.info(f'Waiting to {output_folder} to be generated... (waited {i} seconds)')
        i += 1
        if i > max_waiting_time:
            raise OSError(
                f'{output_folder} was not generated after {max_waiting_time} second. Failed to continue the analysis.')
        sleep(1)


def send_email_in_pipeline_end(logger, process_id, email_address, job_name, state):
    email_addresses = [flask_interface_consts.OWNER_EMAIL]
    email_addresses.extend(flask_interface_consts.ADDITIONAL_OWNER_EMAILS)
    if email_address is not None:
        email_addresses.append(email_address)
    else:
        logger.warning(f'process_id = {process_id} email_address is None, state = {state}, job_name = {job_name}')

    # sends mail once the job finished or crashes
    if state == SharedConsts.State.Finished:
        send_email(logger, SharedConsts.EMAIL_CONSTS.create_title(state, job_name),
                   SharedConsts.EMAIL_CONSTS.CONTENT_PROCESS_FINISHED.format(process_id=process_id), email_addresses)
    elif state == SharedConsts.State.Crashed:
        send_email(logger, SharedConsts.EMAIL_CONSTS.create_title(state, job_name),
                   SharedConsts.EMAIL_CONSTS.CONTENT_PROCESS_CRASHED.format(process_id=process_id), email_addresses)


def add_results_to_final_dir(logger, source, final_output_dir):
    source_dir_name = os.path.split(source)[1]
    dest = os.path.join(final_output_dir, consts.OUTPUTS_DIRECTORIES_MAP[source_dir_name])

    try:
        logger.info(f'Copying {source} TO {dest}')
        shutil.copytree(source, dest, dirs_exist_ok=True)
    except Exception:
        logger.exception(f'Failed to copy {source} to {dest}')
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


def get_job_logger(log_file_dir, level=logging.INFO):
    job_name = os.environ.get(consts.JOB_NAME_ENVIRONMENT_VARIABLE, ''.join(random.choice(string.ascii_letters) for _ in range(5)))
    job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, ''.join(random.choice(string.digits) for _ in range(5)))

    logger = logging.getLogger('main')
    file_handler = logging.FileHandler(os.path.join(log_file_dir, f'{job_name}_{job_id}_log.txt'), mode='a')
    formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(level)

    return logger


def get_job_times_logger(log_file_dir, level=logging.INFO):
    job_name = os.environ.get(consts.JOB_NAME_ENVIRONMENT_VARIABLE, ''.join(random.choice(string.ascii_letters) for _ in range(5)))
    job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, ''.join(random.choice(string.digits) for _ in range(5)))

    logger = logging.getLogger('times')
    file_handler = logging.FileHandler(os.path.join(log_file_dir, f'{job_name}_{job_id}_times_log.txt'), mode='a')
    formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(level)

    return logger


def str_to_bool(s):
    return s.lower() in ['true', '1', 't', 'y', 'yes']


def none_or_str(value):
    if value == 'None':
        return None
    return value