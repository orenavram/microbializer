import subprocess
import os
import logging
from time import time, sleep, ctime
from subprocess import run

logger = logging.getLogger('main')  # use logger instead of printing


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
        #raise ValueError(f'\n{"#"*100}\nnum_of_expected_results is {num_of_expected_results}! Something went wrong in the previous step_number...\n{"#"*100}')
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
    outputs_dir = os.path.join(outputs_dir_prefix, dir_name)
    logger.info(f'{ctime()}: Creating {outputs_dir}\n')
    os.makedirs(outputs_dir, exist_ok=True)

    tmp_dir = os.path.join(tmp_dir_prefix, dir_name)
    logger.info(f'{ctime()}: Creating {tmp_dir}\n')
    os.makedirs(tmp_dir, exist_ok=True)

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


def new_submit_pipeline_step(script_path, params_lists, tmp_dir, job_name, queue_name, q_submitter_script_path,
                             new_line_delimiter, verbose=False, required_modules_as_list=None, num_of_cpus=1):

    required_modules_as_str = 'python/python-anaconda3.6.5'
    if required_modules_as_list:
        # don't forget a space after the python module!!
        required_modules_as_str += ' ' + ' '.join(required_modules_as_list)
    cmds_as_str = f'module load {required_modules_as_str}'
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the rows in q_submitter
    cmds_as_str += new_line_delimiter

    example_cmd = ' '.join(['python', script_path, *[str(param) for param in params_lists[0]]] + (['-v'] if verbose else [])) + ';'
    for params in params_lists:
        cmds_as_str += ' '.join(['python', script_path, *[str(param) for param in params]] + (['-v'] if verbose else [])) + ';'
        cmds_as_str += new_line_delimiter

    cmds_as_str += '\t' + job_name + '\n'
    cmds_path = os.path.join(tmp_dir, f'{job_name}.cmds')
    if os.path.exists(cmds_path):
        cmds_path = os.path.join(tmp_dir, f'{job_name}_{time()}.cmds')

    with open(cmds_path, 'w') as f:
        f.write(cmds_as_str)

    # process_str = f'{q_submitter_script_path} {cmds_path} {tmp_dir} -q {queue_name} --cpu {num_of_cpus}'
    process = [q_submitter_script_path, cmds_path, tmp_dir, '-q', queue_name, '--cpu', str(num_of_cpus)]
    logger.info(f'Calling:\n{" ".join(process)}')
    run(process)
    return example_cmd


def submit_batches(script_path, all_cmds_params, tmp_dir, job_prefix, queue_name, num_of_cmds_per_job=1,
                   q_submitter_script_path='/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py',
                   new_line_delimiter='!@#', required_modules_as_list=None, num_of_cpus=1):

    num_of_batches = 0
    example_cmd_from_batch = 'NO COMMANDS WERE FETCHED'

    for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
        current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
        example_cmd_from_batch = new_submit_pipeline_step(script_path, current_batch, tmp_dir,
                                                          f'{num_of_batches}_{job_prefix}', queue_name,
                                                          q_submitter_script_path, new_line_delimiter, verbose=False,
                                                          required_modules_as_list=required_modules_as_list,
                                                          num_of_cpus=num_of_cpus)
        logger.info(f'Example command from current batch:\n{example_cmd_from_batch}')
        num_of_batches += 1

    return num_of_batches, example_cmd_from_batch



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

# def done_if_exist(done_file, file_to_check = './'):
#     if os.path.exists(file_to_check):
#         execute(['touch', done_file])
#     else:
#         raise AssertionError(file_to_check + ' does not exist....')