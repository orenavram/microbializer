import os
import shutil
import subprocess
from time import time, sleep
import re
import logging
from datetime import timedelta
import pandas as pd

from . import consts
from . import cgi_consts
from .email_sender import send_email
from .q_submitter_power import submit_cmds_from_file_to_q


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


def measure_time(total):
    hours = total // 3600
    minutes = (total % 3600) // 60
    seconds = total % 60
    if hours != 0:
        return f'{hours}:{minutes:02}:{seconds:02} hours'
    elif minutes != 0:
        return f'{minutes}:{seconds:02} minutes'
    else:
        return f'{seconds} seconds'


def execute(logger, process, process_is_string=False):
    process_str = process if process_is_string else ' '.join(str(token) for token in process)
    logger.info(f'Calling (process_is_string == {process_is_string}):\n{process_str}')
    subprocess.run(process, shell=process_is_string)


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
                f'\t{measure_time(total_time)} have passed since started waiting ({num_of_expected_results} - {current_num_of_results} = {jobs_left} more files are still missing)')

    end = time()
    total_time_waited = measure_time(int(end - start))
    logger.info(f'Done waiting for: {script_name} (took {total_time_waited}).')

    if not recursive_step:
        walltime_sum, log_files_without_times = get_jobs_cummulative_time(path)
    else:
        walltime_sum, log_files_without_times, sub_steps_cummulative_times = get_recursive_step_cummulative_times(path)
        for step_name, step_time in sub_steps_cummulative_times.items():
            times_logger.info(f'Sub-step {step_name} took cumulatively {step_time} wallclock time')

    times_logger.info(f'Step {script_name} took {total_time_waited}. '
                      f'There were {num_of_expected_results} jobs and '
                      f'cumulatively they took {walltime_sum} wallclock time. ' +
                      (f'Times are not complete since files {log_files_without_times} do not have times records'
                       if log_files_without_times else ''))

    assert not os.path.exists(error_file_path)


def get_job_time_from_log_file(log_file_content, pattern_for_runtime):
    runtime_string_match = pattern_for_runtime.search(log_file_content)
    if not runtime_string_match:
        return None

    runtime_string = runtime_string_match.group(1)

    if '-' in runtime_string:
        days, time = runtime_string.split('-')
    else:
        days = 0
        time = runtime_string
    hours, minutes, seconds = map(float, time.split(':'))
    runtime = timedelta(days=int(days), hours=hours, minutes=minutes, seconds=seconds)

    return runtime


def get_jobs_cummulative_time(path):
    log_files = [file_path for file_path in os.listdir(path) if file_path.endswith('log.txt')]
    pattern_for_walltime = re.compile(f'{consts.JOB_WALL_TIME_KEY}(.+) TimeLimit')

    walltime_sum = timedelta()
    log_files_without_times = []
    for log_file_name in log_files:
        log_file_full_path = os.path.join(path, log_file_name)
        with open(log_file_full_path, 'r') as log_file:
            content = log_file.read()

            walltime = get_job_time_from_log_file(content, pattern_for_walltime)
            if walltime is None:
                log_files_without_times.append(log_file_full_path)
                continue
            walltime_sum += walltime

    return walltime_sum, log_files_without_times


def get_recursive_step_cummulative_times(path):
    times_log_files = [file_path for file_path in os.listdir(path) if file_path.endswith('times_log.txt')]
    pattern = re.compile(f'Step (.+) took (.+) cumulatively they took (.+) wallclock time')

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
                step_name, step_cummulative_time = match.group(1), match.group(3)
                if step_name in sub_steps_cummulative_times:
                    sub_steps_cummulative_times[step_name] += pd.Timedelta(step_cummulative_time)
                else:
                    sub_steps_cummulative_times[step_name] = pd.Timedelta(step_cummulative_time)

    total_duration = sum(sub_steps_cummulative_times.values(), pd.Timedelta(0))
    return total_duration, log_files_without_times, sub_steps_cummulative_times


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


def submit_mini_batch(logger, script_path, mini_batch_parameters_list, logs_dir, error_file_path, queue_name, account_name=None, job_name='',
                      new_line_delimiter='!@#', verbose=False, required_modules_as_list=None, num_of_cpus=1,
                      memory=None, submit_as_a_job=True, done_file_is_needed=True,
                      done_files_script_path=os.path.join(consts.PROJECT_ROOT_DIR, 'pipeline/auxiliaries/file_writer.py'),
                      command_to_run_before_script=None):
    """
    :param script_path:
    :param mini_batch_parameters_list: a list of lists. each sublist corresponds to a single command and contain its parameters
    :param logs_dir:
    :param job_name:
    :param queue_name:
    :param new_line_delimiter: leave it as is
    :param verbose:
    :param required_modules_as_list: a list of strings containing module names that should be loaded before running
    :param num_of_cpus:
    :param done_files_script_path: leave it as is
    :param submit_as_a_job: if False, fetched directly in shell and not as an independent job
    :return: an example command to debug on the shell
    """

    shell_cmds_as_str = ''

    if not consts.USE_CONDA:
        # COMMAND FOR LOADING RELEVANT MODULES
        required_modules_as_str = 'python/python-anaconda3.6.5'
        if required_modules_as_list:
            # don't forget a space after the python module!!
            required_modules_as_str += ' ' + ' '.join(required_modules_as_list)

        shell_cmds_as_str = f'module load {required_modules_as_str}'
        shell_cmds_as_str += new_line_delimiter  # several commands that will be split to different lines
        # (long lines with ";" are bad practice)
    else:
        # shell_cmds_as_str += f'source ~/.bashrc{new_line_delimiter}'
        conda_sh_path = os.path.join(consts.CONDA_INSTALLATION_DIR, 'etc/profile.d/conda.sh')
        shell_cmds_as_str += f'source {conda_sh_path}{new_line_delimiter}'
        shell_cmds_as_str += f'conda activate {consts.CONDA_ENVIRONMENT_DIR}{new_line_delimiter}'
        shell_cmds_as_str += f'export PATH=$CONDA_PREFIX/bin:$PATH{new_line_delimiter}'

    example_shell_cmd = ' '.join(['python', script_path, *[str(param) for param in mini_batch_parameters_list[0]]] + (
        ['-v'] if verbose else [])) + ';'

    # PREPARING RELEVANT COMMANDS
    if command_to_run_before_script:
        shell_cmds_as_str += f'{command_to_run_before_script}{new_line_delimiter}'

    for params in mini_batch_parameters_list:
        shell_cmds_as_str += ' '.join(
            ['python', script_path, *[str(param) for param in params],
            '-v' if verbose else '', f'--logs_dir {logs_dir}', f'--error_file_path {error_file_path}']) + ';'
        shell_cmds_as_str += new_line_delimiter

    if not job_name:
        job_name = time()

    if done_file_is_needed:
        # GENERATE DONE FILE
        params = [os.path.join(logs_dir, job_name + '.done'), '', f'--logs_dir {logs_dir}']  # write an empty string (like "touch" command)
        shell_cmds_as_str += ' '.join(['python', done_files_script_path, *params]) + ';'
        shell_cmds_as_str += new_line_delimiter

    if submit_as_a_job:
        # WRITING CMDS FILE
        cmds_path = os.path.join(logs_dir, f'{job_name}.cmds')
        with open(cmds_path, 'w') as f:
            f.write(f'{shell_cmds_as_str}\t{job_name}\n')  # ADDING THE JOB NAME

        submit_cmds_from_file_to_q(logger, cmds_path, logs_dir, queue_name, str(num_of_cpus), account_name, new_line_delimiter, memory)
    else:
        # fetch directly on shell
        for shell_cmd in shell_cmds_as_str.split(new_line_delimiter):
            execute(logger, shell_cmd, process_is_string=True)

    return example_shell_cmd


def submit_batch(logger, script_path, batch_parameters_list, logs_dir, error_file_path, job_name_suffix='', queue_name='pupkolabr', account_name=None,
                 num_of_cmds_per_job=1, new_line_delimiter='!@#', required_modules_as_list=None, num_of_cpus=1,
                 memory=None):
    """
    :param script_path:
    :param batch_parameters_list: a list of lists. each sublist corresponds to a single command and contain its parameters
    :param logs_dir:
    :param job_name_suffix: a string that will be concatenated after the batch number as the job name
    :param queue_name:
    :param num_of_cmds_per_job:
    :param new_line_delimiter: leave it as is
    :param required_modules_as_list: a list of strings containing module names that should be loaded before running
    :param num_of_cpus:
    :return: number of batches submitted (in case waiting for the results) and an example command to debug on the shell
    """
    num_of_mini_batches = 0
    example_cmd_from_last_mini_batch = 'NO COMMANDS WERE FETCHED'

    if not job_name_suffix:
        job_name_suffix = time()

    job_name_suffix = job_name_suffix.replace(' ', '_')  # job name cannot contain spaces!

    for i in range(0, len(batch_parameters_list), num_of_cmds_per_job):
        mini_batch_parameters_list = batch_parameters_list[i: i + num_of_cmds_per_job]
        mini_batch_job_name = f'{num_of_mini_batches}_{job_name_suffix}'
        example_cmd_from_last_mini_batch = submit_mini_batch(logger, script_path, mini_batch_parameters_list, logs_dir,
                                                             error_file_path, queue_name, account_name, mini_batch_job_name,
                                                             new_line_delimiter, verbose=False, num_of_cpus=num_of_cpus,
                                                             required_modules_as_list=required_modules_as_list, memory=memory)
        logger.info(f'Example command from batch {mini_batch_job_name}:\n{example_cmd_from_last_mini_batch}')
        num_of_mini_batches += 1

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


def notify_admin(meta_output_dir, meta_output_url, run_number):
    email = 'NO_EMAIL'
    user_email_path = os.path.join(meta_output_dir, cgi_consts.EMAIL_FILE_NAME)
    if os.path.exists(user_email_path):
        with open(user_email_path) as f:
            email = f.read().rstrip()
    error_log_path = 'NO_ERROR_LOG'
    file_with_job_id_on_qstat = os.path.join(meta_output_dir, 'qsub.log')
    if os.path.exists(file_with_job_id_on_qstat):
        with open(file_with_job_id_on_qstat) as f:
            job_id_on_qstat = f.read().strip()
        error_log_path = os.path.join(meta_output_dir, f'{job_id_on_qstat}.ER')
        # TODO: change to ER url and add reading permissions
    # Send me a notification email every time there's a failure
    send_email(smtp_server=consts.SMTP_SERVER,
               sender=consts.ADMIN_EMAIL,
               receiver=consts.OWNER_EMAIL,
               subject=f'{cgi_consts.WEBSERVER_NAME} job {run_number} by {email} has been failed: ',
               content=f"{email}\n\n{os.path.join(meta_output_url, cgi_consts.RESULT_WEBPAGE_NAME)}\n\n"
                       f"{os.path.join(meta_output_url, cgi_consts.CGI_DEBUG_FILE_NAME)}\n\n"
                       f"{os.path.join(meta_output_url, error_log_path)}\n\n"
                       f"{os.path.join(meta_output_dir, error_log_path.replace('ER', 'OU'))}")


def add_results_to_final_dir(logger, source, final_output_dir, keep_in_source_dir=True):
    source_dir_name = os.path.split(source)[1]
    dest = os.path.join(final_output_dir, consts.OUTPUTS_DIRECTORIES_MAP[source_dir_name])

    try:
        if not keep_in_source_dir:
            logger.info(f'Moving {source} TO {dest}')
            shutil.move(source, dest)
        else:
            logger.info(f'Copying {source} TO {dest}')
            shutil.copytree(source, dest)
    except FileExistsError:
        pass

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
    job_name = os.environ.get(consts.JOB_NAME_ENVIRONMENT_VARIABLE, None)
    job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, None)

    logging.basicConfig(format=consts.LOG_MESSAGE_FORMAT,
                        level=level)
    logger = logging.getLogger('main')

    if consts.LOG_IN_SEPARATE_FILES and job_name and job_id:
        file_handler = logging.FileHandler(os.path.join(log_file_dir, f'{job_name}_{job_id}_log.txt'))
        logger.addHandler(file_handler)

    return logger


def get_job_times_logger(log_file_dir, level=logging.INFO):
    job_name = os.environ.get(consts.JOB_NAME_ENVIRONMENT_VARIABLE, None)
    job_id = os.environ.get(consts.JOB_ID_ENVIRONMENT_VARIABLE, None)

    logging.basicConfig(format=consts.LOG_MESSAGE_FORMAT,
                        level=level)
    logger = logging.getLogger('times')

    if consts.LOG_IN_SEPARATE_FILES and job_name and job_id:
        file_handler = logging.FileHandler(os.path.join(log_file_dir, f'{job_name}_{job_id}_times_log.txt'))
        logger.addHandler(file_handler)

    return logger


def str_to_bool(s):
    return s.lower() in ['true', '1', 't', 'y', 'yes']
