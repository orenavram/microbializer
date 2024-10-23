# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 10:16:41 2017

@author: Oren
"""

import os
import time
import subprocess
from . import consts

JOB_EXTENSION = '.slurm'
LOGIN_NODE = 'powerslurm-login'
JOB_SUBMITTER = 'sbatch'
MEMORY_SUFFIX = 'G'
CHECK_JOB_DETAILS_COMMAND = 'scontrol show job'


def add_slurm_header(sbatch_file_content, queue_name, tmp_dir, job_name, CPUs, account_name, memory=None, time_in_hours=None):
    sbatch_file_content += f'#SBATCH --job-name={job_name}\n'
    sbatch_file_content += f'#SBATCH --account={account_name}\n'
    sbatch_file_content += f'#SBATCH --partition={queue_name}\n'
    sbatch_file_content += f'#SBATCH --ntasks=1\n'
    sbatch_file_content += f'#SBATCH --cpus-per-task={CPUs}\n'
    if memory:
        sbatch_file_content += f'#SBATCH --mem={memory}{MEMORY_SUFFIX}\n'
    if time_in_hours:
        sbatch_file_content += f'#SBATCH --time={time_in_hours}:00:00\n'
    sbatch_file_content += f'#SBATCH --output={tmp_dir}/{job_name}_%j.out\n'
    sbatch_file_content += f'#SBATCH --error={tmp_dir}/{job_name}_%j.err\n'

    sbatch_file_content += f'echo Job ID: ${consts.JOB_ID_ENVIRONMENT_VARIABLE}\n'
    sbatch_file_content += f'echo Running on nodes: $SLURM_JOB_NODELIST\n'
    sbatch_file_content += f'echo Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE\n'

    return sbatch_file_content


def generate_job_file(logger, queue_name, tmp_dir, cmds_path, job_name, job_path, CPUs, account_name=None, memory=None,
                      time_in_hours=None):
    """compose the job file content and fetches it"""
    job_file_content = f'#!/bin/bash {"-x" if consts.JOB_FILES_DEBUG_MODE else ""}\n'  # 'old bash: #!/bin/tcsh -x\n'
    job_file_content = add_slurm_header(job_file_content, queue_name, tmp_dir, job_name, CPUs, account_name, memory, time_in_hours)
    job_file_content += f'{cmds_path}\n'

    # log the runtime of the job
    job_log_file_path = f'{tmp_dir}/$(echo ${consts.JOB_NAME_ENVIRONMENT_VARIABLE})_$(echo ${consts.JOB_ID_ENVIRONMENT_VARIABLE})_log.txt'
    job_file_content += f'{CHECK_JOB_DETAILS_COMMAND} ${consts.JOB_ID_ENVIRONMENT_VARIABLE} | grep -m 1 "{consts.JOB_WALL_TIME_KEY}" >> {job_log_file_path}\n'

    with open(job_path, 'w') as job_fp:  # write the job
        job_fp.write(job_file_content)
    subprocess.call(['chmod', '+x', job_path])  # set execution permissions

    logger.debug('First job details for debugging:')
    logger.debug('#' * 80)
    logger.debug('-> job_path is:\n' + job_path)
    logger.debug('\n-> job_file_content is:\n' + job_file_content)
    logger.debug('-> out file is at:\n' + os.path.join(tmp_dir, job_name + '.$JOB_ID.out'))
    logger.debug('#' * 80)


def submit_cmds_from_file_to_q(logger, job_name, cmds_path, tmp_dir, queue_name, CPUs, account_name=None, memory=None,
                               time_in_hours=None, additional_params=''):
    logger.debug('-> Job will be submitted to ' + queue_name + '\'s queue')
    logger.debug('-> out, err and slurm files will be written to:\n' + tmp_dir + '/')
    logger.debug('-> Job will use ' + CPUs + ' CPU(s)\n')

    job_path = os.path.join(tmp_dir, job_name + JOB_EXTENSION)  # path to job
    generate_job_file(logger, queue_name, tmp_dir, cmds_path, job_name, job_path, CPUs, account_name, memory, time_in_hours)

    # execute the job
    # queue_name may contain more arguments, thus the string of the cmd is generated and raw cmd is called

    if consts.Q_SUBMITTER_ADD_SSH_PREFIX:
        terminal_cmd = f'ssh {LOGIN_NODE} "{JOB_SUBMITTER} {job_path} {additional_params}"'  # FIX by danny 5-1-2023
    else:
        terminal_cmd = f'{JOB_SUBMITTER} {job_path} {additional_params}'

    job_submitted_successfully = False
    try_index = 1
    while not job_submitted_successfully:
        try:
            logger.info(f'Submitting: {terminal_cmd}')
            subprocess.run(terminal_cmd, shell=True, capture_output=True, text=True, check=True)
            logger.info(f"Job {job_path} submitted successfully")
            job_submitted_successfully = True
        except subprocess.CalledProcessError as e:
            logger.error(f"Job submission of {job_path} failed (try {try_index}): {e.stderr}")
            try_index += 1
            if try_index >= 100:
                logger.error(f"Job submission of {job_path} failed too many times. Exiting")
                raise e
            time.sleep(1)

    logger.debug('@ -> Sending jobs is done. @')
