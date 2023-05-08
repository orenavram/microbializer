import argparse
import traceback
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def write_to_file(logger, file_path, content=''):  # with a name that is related to the file's name
    if not content:
        content = ''.join(traceback.format_stack())
    with open(file_path, 'w') as f:
        f.write(content)
    logger.info(f'{file_path} was generated.')


if __name__ == '__main__':
    # This block will be executed only when you run it as your main program.
    # If this module is being imported from another script, this block won't be executed, however the function will be available...
    parser = argparse.ArgumentParser()
    parser.add_argument('logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('file_path', help='A path to a file to write to')
    parser.add_argument('--content', help='The content that will be written to the file', default='')
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir)

    write_to_file(logger, args.file_path, args.content)
