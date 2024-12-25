from sys import argv
import argparse
import logging
import subprocess
import os
import sys
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def mcl(logger, input_file, output_file, cpus):
    if os.path.exists(output_file):
        return

    # --abc for a columns format, i.e., item1\item2\tscore
    cmd = f'mcl "{input_file}" -I 1.5 --abc -o "{output_file}" -te {cpus} -V all'
    logger.info(f'Starting MCL. Calling: {cmd}')
    subprocess.run(cmd, shell=True)
    logger.info(f'MCL finished. Output written to {output_file}')


def run_mcl_on_all_putative_ogs(logger, mcl_input_dir, start_og_number, end_og_number, output_dir, cpus):
    for og_number in range(start_og_number, end_og_number + 1):
        input_path = os.path.join(mcl_input_dir, f'OG_{og_number}.mcl_input')
        output_path = os.path.join(output_dir, f'OG_{og_number}.mcl_output')
        mcl(logger, input_path, output_path, cpus)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('mcl_input_dir', help='path to dir of mcl input files')
    parser.add_argument('start_og_number', type=int, help='OG number to start from')
    parser.add_argument('end_og_number', type=int, help='OG number to end at. Inclusive.')
    parser.add_argument('output_dir', help='path to dir the MCL analysis will be written')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        run_mcl_on_all_putative_ogs(logger, args.mcl_input_dir, args.start_og_number, args.end_og_number, args.output_dir, args.cpus)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
