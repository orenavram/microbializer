def prepare(wd_path):

    import os
    error_file = os.path.join(wd_path, 'error.txt')
    try:
        os.remove(error_file)
        print(f'{error_file} was deleted!')
    except:
        print(f'No {error_file} file to delete was found.')

    html_path = os.path.join(wd_path, 'output.html')
    with open(html_path) as f:
        html_content = ''
        for line in f:
            html_content += line
            if line.startswith('<!--result-->'):
                break

    import sys
    sys.path.insert(0, '/bioseq/microbializer/pipeline')
    from CONSTANTS import RELOAD_TAGS
    html_content = html_content.replace('FAILED', 'RUNNING')
    html_content = html_content.replace('FINISHED', 'RUNNING')
    html_content = html_content.replace('progress-bar-striped', 'progress-bar-striped active')
    html_content = html_content.replace(f'<!--{RELOAD_TAGS}-->', RELOAD_TAGS)

    with open(html_path, 'w') as f:
        f.write(html_content)

    print(f'{html_path} was reverted!')


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('wd_path', help="A path to a working dir of microbializer's job.\n"
                                        "E.g., /bioseq/data/results/microbializer/158875031326946667844691750504")
    args = parser.parse_args()

    prepare(args.wd_path)