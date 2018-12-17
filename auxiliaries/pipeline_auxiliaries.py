import logging
logger = logging.getLogger('main') # use logger instead of printing


import subprocess
import os
from auxiliaries.directory_creator import create_dir
from auxiliaries.email_sender import send_email
from time import time, sleep


def measure_time(total):
    hours = total // 3600
    minutes = (total% 3600) // 60
    seconds = total % 60
    if hours != 0:
        return f'{hours}:{minutes}:{seconds} hours'
    elif minutes != 0:
        return f'{minutes}:{seconds} minutes'
    else:
        return f'{seconds} seconds'


def execute(process, raw=False):
    process_str = process if raw else ' '.join(str(token) for token in process)
    logger.warning(f'Calling: {process_str} (raw == {raw})')
    subprocess.call(process, shell=raw)


def wait_for_results(script_name, path, num_of_expected_results, suffix='done', remove=False, time_to_wait=100):
    '''waits until path contains num_of_expected_results $suffix files'''
    start = time()
    logger.warning(f'Waiting for {script_name}... (continues when {num_of_expected_results} results will be in {path})')
    if num_of_expected_results==0:
        raise ValueError(f'\n{"#"*50}\nnum_of_expected_results is {num_of_expected_results}! Something went wrong in the previous step...\n{"#"*50}')
    i = 0
    current_num_of_results = 0
    while num_of_expected_results > current_num_of_results:
        current_num_of_results = sum(1 for x in os.listdir(path) if x.endswith(suffix))
        jobs_left = num_of_expected_results - current_num_of_results
        sleep(time_to_wait)
        i += time_to_wait
        logger.warning(f'{i} seconds have passed since started waiting ({num_of_expected_results} - {current_num_of_results} = {jobs_left} more files are still missing)')
    if remove:
        execute(['python', '-u', '/groups/pupko/orenavr2/pipeline/RemoveDoneFiles.py', path, suffix])
    end = time()
    logger.warning(f'Done waiting for:\n{script_name}\n(took {measure_time(int(end-start))}).\n')


# def remove_files_with_suffix(path, suffix='done'):
#     '''remove all files from path that end with suffix'''
#     logger.warning(f'Removing {suffix} files from {path}')
#     for file_name in os.listdir(path):
#         if file_name.endswith(suffix):
#             file_path = os.path.join(path,file_name)
#             logger.debug(f'Removing {file_path}')
#             os.remove(file_path)
#     logger.warning('Done removing.')


def prepare_directories(outputs_dir_prefix, tmp_dir_prefix, dir_name):
    outputs_dir = os.path.join(outputs_dir_prefix, dir_name)
    create_dir(outputs_dir)
    tmp_dir = os.path.join(tmp_dir_prefix, dir_name)
    create_dir(tmp_dir)
    return outputs_dir, tmp_dir


def submit_pipeline_step(script_path, params, tmp_dir, job_name, queue_name, new_line_delimiter='!@#', q_submitter_script_path='/bioseq/bioSequence_scripts_and_constants/q_submitter.py', done_files_script_path='/groups/pupko/orenavr2/microbializer/auxiliaries/file_writer.py', required_modules = []):

    required_modules_as_str = ' '.join(['python/anaconda_python-3.6.4'] + required_modules)
    cmds = f'module load {required_modules_as_str}'
    cmds += new_line_delimiter

    # ACTUAL COMMAND
    cmds += ' '.join(['python', '-u', script_path, *params, ';'])
    cmds += new_line_delimiter # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the rows in q_submitter

    # GENERATE DONE FILE
    params = [os.path.join(tmp_dir, job_name + '.done'), ''] # write an empty string (like "touch" command)
    cmds += ' '.join(['python', '-u', done_files_script_path, *params, ';'])
    cmds += new_line_delimiter

    cmds += '\t' + job_name
    cmds_path = os.path.join(tmp_dir, job_name + '.cmds')
    with open(cmds_path, 'w') as f:
        f.write(cmds)
    execute([q_submitter_script_path, cmds_path, tmp_dir, '-q', queue_name])

# def done_if_exist(done_file, file_to_check = './'):
#     if os.path.exists(file_to_check):
#         execute(['touch', done_file])
#     else:
#         raise AssertionError(file_to_check + ' does not exist....')