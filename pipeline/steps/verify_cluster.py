from sys import argv
import argparse
import logging
import os
import sys
import shutil
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def verify(logger, input_file, output_dir):
    input_og_name = os.path.splitext(os.path.basename(input_file))[0]
    with open(input_file, 'r') as f:
        lines = [line.rstrip() for line in f if line]

    if len(lines) == 0:
        raise ValueError(f'{input_file} is empty! There\'s a bug in the previous step!')
    elif len(lines) == 1:
        output_file_path = os.path.join(output_dir, input_og_name + ".verified_cluster")
        if os.path.exists(output_file_path):
            return
        shutil.copyfile(input_file, output_file_path)
    else:  # 1 < len(lines)
        og_subset_id = 0
        for line in lines:
            if len(line.split('\t')) == 1:
                continue  # Skip lines with 1 protein (orphan genes)
            verified_cluster_path = os.path.join(output_dir, f"{input_og_name}_{og_subset_id}.split_cluster")
            if not os.path.exists(verified_cluster_path):
                with open(verified_cluster_path, 'w') as verified_cluster_file:
                    verified_cluster_file.write(line)
            og_subset_id += 1


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='path to an MCL analysis file')
    parser.add_argument('output_dir',
                        help='dir path to which the MCL analysis will be moved if clustering criterion was met')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        verify(logger, args.input_file, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
