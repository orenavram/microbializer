from pathlib import Path
import argparse
import shutil
import os
import logging
import pandas as pd
import math

from . import consts


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ("yes", "true", "t", "1"):
        return True
    elif value.lower() in ("no", "false", "f", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def none_or_path(value):
    if value == 'None':
        return None
    return Path(value)


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


def fail(logger, error_msg, error_file_path):
    logger.error(error_msg)
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise ValueError(error_msg)


def write_done_file(logger, done_file_path):
    with open(done_file_path, 'w') as f:
        f.write('.')
    logger.info(f'{done_file_path} was generated.')


def get_logger(log_file_path, logger_name, verbose: bool):
    logger = logging.getLogger(logger_name)
    file_handler = logging.FileHandler(log_file_path, mode='a')
    formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    return logger


def flatten(l):
    return [item.strip() for sublist in l for item in sublist if pd.notna(item)]


def get_directory_size_in_gb(directory):
    total_size = 0
    for file_path in directory.iterdir():
        # Skip broken symlinks and non-files
        if file_path.is_file():
            total_size += file_path.stat().st_size
    # Convert bytes to gigabytes and round up
    size_in_gb = total_size / (1024 ** 3)
    return size_in_gb


def get_required_memory_gb_to_load_csv(csv_path: Path):
    csv_size_bytes = csv_path.stat().st_size
    requited_memory_bytes = csv_size_bytes * 10
    requited_memory_gb = math.ceil(requited_memory_bytes / (1024 ** 3))
    return requited_memory_gb
