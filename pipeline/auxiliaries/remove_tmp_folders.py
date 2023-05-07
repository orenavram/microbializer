def remove_dirs_from_tmp_dir(tmp_dir):
    import shutil
    import os
    for file in os.listdir(tmp_dir):
        folder = os.path.join(tmp_dir, file)
        if os.path.isdir(folder):
            logger.info(f'Removing {folder}')
            shutil.rmtree(folder, ignore_errors=True)


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('tmp_dir', help='A path from which all folders will be removed')
    args = parser.parse_args()

    import logging

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    remove_dirs_from_tmp_dir(args.tmp_dir)
