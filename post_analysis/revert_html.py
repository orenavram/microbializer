def revert_html(html_path):
    with open(html_path) as f:
        html_content = ''
        for line in f:
            html_content += line
            if line.startswith('<!--result-->'):
                break

    from sys import path
    path.insert(0, '/bioseq/microbializer/pipeline')
    from CONSTANTS import RELOAD_TAGS
    html_content = html_content.replace('FAILED', 'RUNNING')
    html_content = html_content.replace('FINISHED', 'RUNNING')
    html_content = html_content.replace('progress-bar-striped', 'progress-bar-striped active')
    html_content = html_content.replace(f'<!--{RELOAD_TAGS}-->', RELOAD_TAGS)

    with open(html_path, 'w') as f:
        f.write(html_content)

if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('html_path', help='path an output.html file')
    args = parser.parse_args()

    revert_html(args.html_path)