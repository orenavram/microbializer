import subprocess
import os
from directory_creator import create_dir
from time import time, sleep, ctime
import logging

logger = logging.getLogger('main') # use logger instead of printing

def measure_time(total):
    hours = total // 3600
    minutes = (total% 3600) // 60
    seconds = total % 60
    if hours != 0:
        return f'{hours}:{minutes:02}:{seconds:02} hours'
    elif minutes != 0:
        return f'{minutes}:{seconds:02} minutes'
    else:
        return f'{seconds} seconds'


def execute(process, raw=False):
    process_str = process if raw else ' '.join(str(token) for token in process)
    logger.info(f'Calling (raw == {raw}):\n{process_str}')
    subprocess.run(process, shell=raw)


def wait_for_results(script_name, path, num_of_expected_results, error_file_path, suffix='done',
                     remove=False, time_to_wait=10, start=0):
    '''waits until path contains num_of_expected_results $suffix files'''
    if not start:
        start = time()
    logger.info(f'Waiting for {script_name}...\nContinues when {num_of_expected_results} results will be in:\n{path}')
    if num_of_expected_results==0:
        logger.fatal(f'\n{"#"*100}\nnum_of_expected_results in {path} is {num_of_expected_results}!\nSomething went wrong in the previous step...\n{"#"*100}')
        #raise ValueError(f'\n{"#"*100}\nnum_of_expected_results is {num_of_expected_results}! Something went wrong in the previous step...\n{"#"*100}')
    total_time = 0
    i = 0
    current_num_of_results = 0
    while num_of_expected_results > current_num_of_results:
        current_num_of_results = sum(1 for x in os.listdir(path) if x.endswith(suffix))
        jobs_left = num_of_expected_results - current_num_of_results
        sleep(time_to_wait)
        total_time += time_to_wait
        i += 1
        if i % 5 == 0:  # print status every 5 cycles of $time_to_wait
            logger.info(f'\t{measure_time(total_time)} have passed since started waiting ({num_of_expected_results} - {current_num_of_results} = {jobs_left} more files are still missing)')
    # if remove:
    #     execute(['python', '-u', '/groups/pupko/orenavr2/pipeline/RemoveDoneFiles.py', path, suffix])
    end = time()
    logger.info(f'Done waiting for:\n{script_name}\n(took {measure_time(int(end-start))}).\n')
    assert not os.path.exists(error_file_path)


# def remove_files_with_suffix(path, suffix='done'):
#     '''remove all files from path that end with suffix'''
#     logger.info(f'Removing {suffix} files from {path}')
#     for file_name in os.listdir(path):
#         if file_name.endswith(suffix):
#             file_path = os.path.join(path,file_name)
#             logger.debug(f'Removing {file_path}')
#             os.remove(file_path)
#     logger.info('Done removing.')


def prepare_directories(outputs_dir_prefix, tmp_dir_prefix, dir_name):
    logger.info(ctime())
    outputs_dir = os.path.join(outputs_dir_prefix, dir_name)
    create_dir(outputs_dir)
    tmp_dir = os.path.join(tmp_dir_prefix, dir_name)
    create_dir(tmp_dir)
    return outputs_dir, tmp_dir


def submit_pipeline_step(script_path, params, tmp_dir, job_name, queue_name, new_line_delimiter='!@#',
                         q_submitter_script_path='/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py',
                         done_files_script_path='/bioseq/microbializer/auxiliaries/file_writer.py',
                         required_modules_as_list=None, more_cmds=None, num_of_cpus=1):

    for char in ' ,;()"\'':
        job_name = job_name.replace(char, '_')

    required_modules_as_str = 'python/python-anaconda3.6.5-orenavr2'
    if required_modules_as_list:
        # don't forget a space after the python module!!
        required_modules_as_str += ' ' + ' '.join(required_modules_as_list)
    cmds_as_str = f'module load {required_modules_as_str}'
    cmds_as_str += new_line_delimiter

    if more_cmds:
        for cmd in more_cmds:
            #cmds_as_str += ' '.join(['python', '-u', script_path, *cmd]) + ';' #unbuffering
            cmds_as_str += ' '.join(['python', script_path, *cmd]) + ';' #buffering
            cmds_as_str += new_line_delimiter  # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the rows in q_submitter

    # ACTUAL COMMAND (last command if it's a batch)
    #cmds_as_str += ' '.join(['python', '-u', script_path, *params])+';' #unbuffering
    cmds_as_str += ' '.join(['python', script_path, *params])+';' #buffering
    cmds_as_str += new_line_delimiter # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the rows in q_submitter

    # GENERATE DONE FILE
    params = [os.path.join(tmp_dir, job_name + '.done'), ''] # write an empty string (like "touch" command)
    cmds_as_str += ' '.join(['python', done_files_script_path, *params])+';'
    cmds_as_str += new_line_delimiter

    cmds_as_str += '\t' + job_name + '\n'
    logger.debug(cmds_as_str)
    cmds_path = os.path.join(tmp_dir, job_name + '.cmds')
    with open(cmds_path, 'w') as f:
        f.write(cmds_as_str)
    execute([q_submitter_script_path, cmds_path, tmp_dir, '-q', queue_name, '--cpu', str(num_of_cpus)])


def fail(error_msg, error_file_path):
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise ValueError(error_msg)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

# def done_if_exist(done_file, file_to_check = './'):
#     if os.path.exists(file_to_check):
#         execute(['touch', done_file])
#     else:
#         raise AssertionError(file_to_check + ' does not exist....')