# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 10:16:41 2017

@author: Oren
"""

import argparse
import logging
import os
from subprocess import call
from . import consts


def generate_qsub_file(logger, queue_name, tmp_dir, cmd, prefix_name, qsub_path, CPUs):
    """compose qsub_file content and fetches it"""
    qsub_file_content = '#!/bin/bash -x\n'  # 'old bash: #!/bin/tcsh -x\n'
    qsub_file_content += '#PBS -S /bin/bash\n'  # '#PBS -S /bin/tcsh\n'
    if int(CPUs) > 1:
        qsub_file_content += f'#PBS -l ncpus={CPUs}\n'
    qsub_file_content += f'#PBS -q {queue_name}@power9\n'
    qsub_file_content += f'#PBS -N {prefix_name}\n'
    qsub_file_content += f'#PBS -e {tmp_dir}\n'  # error log
    qsub_file_content += f'#PBS -o {tmp_dir}\n'  # output log
    qsub_file_content += f'#PBS -r y\n'  # makes the job Rerunnable
    qsub_file_content += f'hostname\n'
    qsub_file_content += f'echo job_name: {prefix_name}\n'
    qsub_file_content += f'echo $PBS_JOBID\n'
    qsub_file_content += f'{cmd}\n'
    with open(qsub_path, 'w') as f_qsub:  # write the job
        f_qsub.write(qsub_file_content)
    call(['chmod', '+x', qsub_path])  # set execution permissions

    logger.debug('First job details for debugging:')
    logger.debug('#' * 80)
    logger.debug('-> qsub_path is:\n' + qsub_path)
    logger.debug('\n-> qsub_file_content is:\n' + qsub_file_content)
    logger.debug('-> out file is at:\n' + os.path.join(tmp_dir, prefix_name + '.$JOB_ID.out'))
    logger.debug('#' * 80)


def submit_cmds_from_file_to_q(logger, cmds_file, tmp_dir, queue_name, CPUs, dummy_delimiter='!@#', start=0, end=float('inf'),
                               additional_params=''):
    logger.debug('-> Jobs will be submitted to ' + queue_name + '\'s queue')
    logger.debug('-> out, err and pbs files will be written to:\n' + tmp_dir + '/')
    logger.debug('-> Jobs are based on cmds ' + str(start) + ' to ' + str(end) + ' (excluding) from:\n' + cmds_file)
    logger.debug('-> Each job will use ' + CPUs + ' CPU(s)\n')

    logger.debug('Starting to send jobs...')
    cmd_number = 0
    with open(cmds_file) as f_cmds:
        for line in f_cmds:
            if int(start) <= cmd_number < end and not line.isspace():
                try:
                    cmd, prefix_name = line.rstrip().split('\t')
                except:
                    logger.error(f'UNABLE TO PARSE LINE:\n{line}')
                    raise
                    # the queue does not like very long commands so I use a dummy delimiter (!@# by default) to break the rows:
                cmd = cmd.replace(dummy_delimiter, '\n')
                qsub_path = os.path.join(tmp_dir, prefix_name + '.pbs')  # path to job

                generate_qsub_file(logger, queue_name, tmp_dir, cmd, prefix_name, qsub_path, CPUs)

                # execute the job
                # queue_name may contain more arguments, thus the string of the cmd is generated and raw cmd is called

                if consts.TEST:
                    terminal_cmd = f'/opt/pbs/bin/qsub {qsub_path} {additional_params}'
                else:
                    terminal_cmd = f'ssh power9login "/opt/pbs/bin/qsub {qsub_path} {additional_params}"'  # FIX by danny 5-1-2023
                logger.info(f'Submitting: {terminal_cmd}')

                call(terminal_cmd, shell=True)

            cmd_number += 1

    logger.debug('@ -> Sending jobs is done. @')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('cmds_file',
                        help='A file containing jobs commands to execute in the queue. Each row should contain a (set of) command(s separated by $dummy_delimiter) then a "\t" and a job name',
                        type=lambda file_path: str(file_path) if os.path.exists(file_path) else parser.error(
                            f'{file_path} does not exist!'))
    parser.add_argument('tmp_dir', help='A temporary directory where the log files will be written to')
    parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                        default='pupkolab')  # , choices=['pupko', 'itaym', 'lilach', 'bioseq'])
    parser.add_argument('--cpu', help='How many CPUs will be used?', choices=[str(i) for i in range(1, 29)],
                        default='1')
    parser.add_argument('--dummy_delimiter',
                        help='The queue does not "like" very long commands; A dummy delimiter is used to break each row into different commands of a single job',
                        default='!@#')
    parser.add_argument('--start', help='Skip jobs until $start', type=int, default=0)
    parser.add_argument('--end', help='Skip jobs from $end+1', type=int, default=float('inf'))
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--additional_params', help='Other specific parameters, such as, which machine to use',
                        default='')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logger = logging.getLogger('main')
    logger.debug(f'args = {args}')

    if not os.path.exists(args.tmp_dir):
        logger.debug(f'{args.tmp_dir} does not exist. Creating tmp path...')
        os.makedirs(args.tmp_dir, exist_ok=True)

    submit_cmds_from_file_to_q(logger, args.cmds_file, args.tmp_dir, args.queue_name, args.cpu, args.dummy_delimiter,
                               args.start, args.end, args.additional_params)
    # print(args.cmds_file, args.tmp_dir, args.queue_name, args.cpu, args.verbose, args.start, args.end)
