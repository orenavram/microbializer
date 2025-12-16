import os
import subprocess
from time import time, sleep
import re
from datetime import timedelta
import pandas as pd
import stat
from pathlib import Path
import traceback
from sys import argv
from concurrent.futures import ProcessPoolExecutor

from . import consts
from .q_submitter_power import submit_cmds_from_file_to_q
from .general_utils import get_logger, fail, str_to_bool


def validate_slurm_error_logs(logger, slurm_logs_dir, error_file_path):
    for slurm_log_file in slurm_logs_dir.glob('*.err'):
        if slurm_log_file.stat().st_size == 0:
            continue

        with open(slurm_log_file) as f:
            for line in f:
                if 'error: Detected 1 oom_kill event' in line:  # Catch memory errors
                    fail(logger, f'file {slurm_log_file} shows a slurm error: {line}', error_file_path)


def wait_for_results(logger, times_logger, script_name, path, error_file_path, start=0, recursive_step=False):
    """waits until path contains num_of_expected_results .done files"""
    if (path / consts.STEP_INPUTS_DIR_NAME).exists():
        num_of_expected_results = sum(1 for _ in (path / consts.STEP_INPUTS_DIR_NAME).glob('*.txt'))
    else:
        num_of_expected_results = 1

    logger.info(f'Waiting for {script_name}... Continues when {num_of_expected_results} results will be in: {path}')

    if num_of_expected_results == 0:
        raise ValueError('Number of expected results is 0! Something went wrong in the previous analysis steps.')

    if not start:
        start = time()

    i = 0
    current_num_of_results = 0
    while num_of_expected_results > current_num_of_results:
        assert not error_file_path.exists()

        current_num_of_results = sum(1 for _ in path.glob(f'*{consts.JOB_DONE_FILE_SUFFIX}'))
        jobs_left = num_of_expected_results - current_num_of_results
        sleep(consts.CHECK_JOB_DONE_INTERVAL_SECONDS)
        total_time_waited = timedelta(seconds=int(time() - start))
        i += 1
        if i % 5 == 0:  # print status every 5 cycles of $time_to_wait
            logger.info(
                f'\t{total_time_waited} have passed since started waiting ('
                f'{num_of_expected_results} - {current_num_of_results} = {jobs_left} more files are still missing)')

    total_time_waited = timedelta(seconds=int(time() - start))
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
    elif 'day' in runtime_string:
        days_string, time_string = runtime_string.split(' day, ')
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


def submit_job(logger, config, script_path, script_parameters, logs_dir, job_name, num_of_cpus=1,
               memory=None, time_in_hours=None, environment_variables_to_change_before_script: dict = None,
               alternative_error_file=None, job_array_interval=None, job_input_path=None):
    """
    :param script_path:
    :param script_parameters: a list of parameters
    :param logs_dir:
    :param job_name:
    :param num_of_cpus:
    :return: an example command to debug on the shell
    """

    shell_cmds_as_str = ''

    if config.use_job_manager:
        # shell_cmds_as_str += f'source ~/.bashrc{new_line_delimiter}'
        conda_sh_path = config.conda_installation_dir / 'etc' / 'profile.d' / 'conda.sh'
        shell_cmds_as_str += f'source {conda_sh_path}\n'
        shell_cmds_as_str += f'conda activate {config.conda_environment_dir}\n'
        shell_cmds_as_str += f'export PATH=$CONDA_PREFIX/bin:$PATH\n'

    # PREPARING RELEVANT COMMANDS
    if environment_variables_to_change_before_script and config.use_job_manager:
        for key, value in environment_variables_to_change_before_script.items():
            shell_cmds_as_str += f'export {key}={value}\n'

    error_file_path = alternative_error_file if alternative_error_file else config.error_file_path

    shell_cmds_as_str += ' '.join(
        ['python', str(script_path), *[str(param) for param in script_parameters],
         f'-v {config.verbose}', f'--logs_dir {logs_dir}', f'--error_file_path {error_file_path}',
         f'--job_name {job_name}', f'--use_job_manager {config.use_job_manager}', f'--cpus {num_of_cpus}'])

    if job_input_path:
        shell_cmds_as_str += f' --job_input_path {job_input_path}'

    # GENERATE DONE FILE
    done_file_name = f'{job_name}_${consts.JOB_ID_ENVIRONMENT_VARIABLE}{consts.JOB_DONE_FILE_SUFFIX}'
    shell_cmds_as_str += f'\ntouch {logs_dir / done_file_name}\n'

    # WRITING CMDS FILE
    cmds_path = logs_dir / f'{job_name}.sh'
    with open(cmds_path, 'w') as f:
        f.write(shell_cmds_as_str)

    if config.use_job_manager:
        # Add execution permissions to cmds_path
        current_permissions = cmds_path.stat().st_mode
        cmds_path.chmod(current_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        if memory is None:
            memory = config.job_default_memory_gb

        submit_cmds_from_file_to_q(logger, job_name, cmds_path, logs_dir, config.queue_name, config.qos, str(num_of_cpus),
                                   config.account_name, memory, time_in_hours, config.node_name, job_array_interval)
    else:
        new_env = os.environ.copy()
        if environment_variables_to_change_before_script:
            for key, value in environment_variables_to_change_before_script.items():
                new_env[key] = value
                logger.info(f'Changing environment variable {key} to {value} in the environment of {script_path}')

        # fetch directly on shell
        for shell_cmd in shell_cmds_as_str.split('\n'):
            if shell_cmd:
                logger.info(f'Running command: {shell_cmd}')
                subprocess.run(shell_cmd, shell=True, check=True, capture_output=True, text=True, env=new_env)


def submit_batch(logger, config, script_path, script_parameters, logs_dir, job_name_suffix,
                 num_of_cpus=1, memory=None, time_in_hours=None):
    """
    :param script_path:
    :param batch_parameters_list: a list of lists. each sublist corresponds to a single command and contain its parameters
    :param logs_dir:
    :param job_name_suffix: a string that will be concatenated after the batch number as the job name
    :param num_of_cpus:
    :return: number of batches submitted (in case waiting for the results) and an example command to debug on the shell
    """
    job_name_suffix = job_name_suffix.replace(' ', '_')  # job name cannot contain spaces!
    batch_inputs_dir = logs_dir / consts.STEP_INPUTS_DIR_NAME

    if (config.use_job_manager and not config.use_job_array) or config.max_parallel_jobs == 1:
        for i, input_file in enumerate(batch_inputs_dir.glob('*.txt')):
            job_name = f'{i}_{job_name_suffix}'
            submit_job(
                logger, config, script_path, script_parameters, logs_dir, job_name,
                num_of_cpus=num_of_cpus, memory=memory, time_in_hours=time_in_hours, job_input_path=input_file)
            sleep(0.1)

    elif config.use_job_manager and config.use_job_array:
        job_count = sum(1 for _ in batch_inputs_dir.glob("*.txt"))
        submit_job(
            logger, config, script_path, script_parameters, logs_dir, job_name_suffix,
            num_of_cpus=num_of_cpus, memory=memory, time_in_hours=time_in_hours,
            job_array_interval=f'0-{job_count - 1}',
            job_input_path=f'{batch_inputs_dir}/${consts.JOB_ARRAY_TASK_ID_ENVIRONMENT_VARIABLE}.txt')
        sleep(0.1)

    else:
        executor = ProcessPoolExecutor(max_workers=config.max_parallel_jobs)
        for i, input_file in enumerate(batch_inputs_dir.glob('*.txt')):
            job_name = f'{i}_{job_name_suffix}'
            executor.submit(submit_job, logger, config, script_path, script_parameters, logs_dir,
                            job_name, num_of_cpus=num_of_cpus, memory=memory, time_in_hours=time_in_hours, job_input_path=input_file)
            sleep(0.1)


def get_job_logger(log_file_dir, job_name, verbose, use_job_manager):
    if use_job_manager:
        job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, '')
    else:
        job_id = ''

    logger = get_logger(log_file_dir / f'{job_name}_{job_id}_log.txt', 'main', verbose)

    return logger


def get_job_times_logger(log_file_dir, job_name, verbose, use_job_manager):
    if use_job_manager:
        job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, '')
    else:
        job_id = ''

    logger = get_logger(log_file_dir / f'{job_name}_{job_id}_times_log.txt', 'times', verbose)

    return logger


def add_default_step_args(args_parser):
    args_parser.add_argument('-v', '--verbose', help='Increase output verbosity', type=str_to_bool)
    args_parser.add_argument('--logs_dir', type=Path, help='path to tmp dir to write logs to')
    args_parser.add_argument('--error_file_path', type=Path, help='path to error file')
    args_parser.add_argument('--job_name', help='job name')
    args_parser.add_argument('--use_job_manager', help='use job manager', type=str_to_bool)
    args_parser.add_argument('--cpus', default=1, type=int)
    args_parser.add_argument('--job_input_path', help='path to an input file for the job', type=Path, required=False)


def run_step(args, step_method, *step_args):
    start_time = time()

    logger = get_job_logger(args.logs_dir, args.job_name, args.verbose, args.use_job_manager)
    logger.info(f'Starting command is: {" ".join(argv)}')

    try:
        step_method(logger, *step_args)
    except subprocess.CalledProcessError as e:
        error_message = f'Error in function "{step_method.__name__}" in command: "{e.cmd}": {e.stderr}'
        logger.exception(error_message)
        with open(args.error_file_path, 'a+') as f:
            f.write(error_message)
    except Exception as e:
        logger.exception(f'Error in function "{step_method.__name__}"')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)

    total_time = timedelta(seconds=int(time() - start_time))
    
    if not args.use_job_manager:
        logger.info(f'Finished running step {step_method.__name__}: {consts.JOB_WALL_TIME_KEY}{total_time} TimeLimit')
        logger.info(f'Finished running step {step_method.__name__}: {consts.JOB_CPUS_KEY}{args.cpus} NumTasks')
