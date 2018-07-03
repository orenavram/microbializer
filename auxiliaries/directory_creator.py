import logging
logger = logging.getLogger('main') # use logger instead of printing

def create_dir(path):
    import os
    if not os.path.exists(path):
        logger.info(f'Creating directory: {path}')
        os.makedirs(path)
    else:
        logger.info(f'Directory already exists: {path}')


if __name__ == '__main__':
    # This block will be executed only when you run it as your main program.
    # If this module is being imported from another script, this block won't be executed, however the function will be available...
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='A path to create')
    args = parser.parse_args()
    create_dir(args.path)