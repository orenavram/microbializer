import logging

logger = logging.getLogger('main')  # use logger instead of printing


def write_to_file(file_path, content=''):  # with a name that is related to the file's name
    if not content:
        import traceback
        content = ''.join(traceback.format_stack())
    with open(file_path, 'w') as f:
        f.write(content)
    logger.debug(f'{file_path} was generated.')


if __name__ == '__main__':
    # This block will be executed only when you run it as your main program.
    # If this module is being imported from another script, this block won't be executed, however the function will be available...
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('file_path', help='A path to a file to write to')
    parser.add_argument('--content', help='The content that will be written to the file', default='')
    args = parser.parse_args()
    write_to_file(args.file_path, args.content)
