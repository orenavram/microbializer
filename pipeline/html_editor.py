import logging
import os
from auxiliaries import consts

logger = logging.getLogger('main')


def add_closing_html_tags(html_path, run_number):
    with open(html_path, 'a') as f:
        f.write(
            f'<hr>\n<h4 class=footer><p align=\'center\'>Questions and comments are welcome! Please '
            f'<span class="admin_link">'
            f'<a href="mailto:{consts.ADMIN_EMAIL}?subject={consts.WEBSERVER_NAME}%20Run%20Number%20{run_number}">contact us</a>'
            f'</span></p></h4>\n'
            f'<div id="bottom_links" align="center"><span class="bottom_link">'
            f'<a href="{consts.WEBSERVER_URL}/" target="_blank">Home</a>'
            f'&nbsp;|&nbsp<a href="{consts.WEBSERVER_URL}/overview.html" target="_blank">Overview</a>\n'
            f'</span>\n'
            f'<br><br><br>\n</body>\n</html>\n')
        f.flush()

    # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
    from time import sleep
    sleep(2 * consts.RELOAD_INTERVAL)
    with open(html_path) as f:
        html_content = f.read()
    html_content = html_content.replace(consts.RELOAD_TAGS, f'<!--{consts.RELOAD_TAGS}-->')
    with open(html_path, 'w') as f:
        f.write(html_content)


def get_html_string_of_result(final_output_dir_name, meta_output_dir, end_of_str, figure_str_to_show_on_html='',
                              raw_str_to_show_on_html='raw data', additional_text=''):
    result = '<tr><td>'
    raw_file_suffix = os.path.join(final_output_dir_name, end_of_str)
    if os.path.exists(os.path.join(meta_output_dir, raw_file_suffix)):
        if figure_str_to_show_on_html:
            result += f'<a href="{raw_file_suffix.replace("txt", "png")}" target="_blank">{figure_str_to_show_on_html}</a> ; ('
        result += f'<a href="{raw_file_suffix}" target="_blank">{raw_str_to_show_on_html}</a>'
        if figure_str_to_show_on_html:
            result += ')\n'
        if additional_text:
            result += additional_text + '\n'
        result += f'<br></td></tr>'

    return result


def edit_success_html(html_path, meta_output_dir, final_output_dir_name, run_number):
    if consts.IGNORE_HTML:
        return

    html_text = ''
    try:
        with open(html_path) as f:
            html_text = f.read()
        # The initial file exists (generate by the cgi) so we can read and parse it.
        html_text = html_text.replace('RUNNING', 'FINISHED').replace(
            f'{consts.WEBSERVER_NAME} is now processing your request. This page will be automatically updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. ',
            '')
    except FileNotFoundError:
        logger.warning(f"Couldn't find html prefix at: {html_path}")

    html_text += f'<div class="container" style="{consts.CONTAINER_STYLE}">\n' \
                 f'<h2>RESULTS:<h2>' \
                 f'<h3><b><a href=\'{final_output_dir_name}.zip\' target=\'_blank\'>Download zipped full results (textual & visual)</a></b></h3>' \
                 f'<table class="table">\n' \
                 f'<thead>\n' \
                 f'<tr><th><h3>Quick access to selected results:</h3></th></tr>\n' \
                 f'</thead>\n' \
                 f'<tbody>'

    raw_file_suffix = os.path.join(final_output_dir_name, '16_species_phylogeny/final_species_tree.txt')
    if os.path.exists(os.path.join(meta_output_dir, raw_file_suffix)):
        html_text += f'<tr><td><a href="{consts.WEBSERVER_URL}/PhyD3/view_tree.php?id={run_number}&f=newick" target="_blank">Interactive species tree</a> ;' \
                     f' (<a href="{raw_file_suffix}" target="_blank">raw data</a>)\n<br></td></tr>'
    else:
        html_text += f'<tr><td>' \
                     f'Species tree is not available for current analysis ' \
                     f'(<a href="https://microbializer.tau.ac.il/faq.html#no_tree" target="_blank">Why?</a>)\n' \
                     f'<br></td></tr>'

    html_text += get_html_string_of_result(final_output_dir_name,
                                           meta_output_dir,
                                            '19_groups_sizes_frequency/groups_sizes_frequency.txt',
                                           figure_str_to_show_on_html='Orthologs groups size dispersion')

    html_text += get_html_string_of_result(final_output_dir_name,
                                           meta_output_dir,
                                            '20_orfs_plots/orfs_counts.txt',
                                           figure_str_to_show_on_html='ORFs per genome dispersion')

    html_text += get_html_string_of_result(final_output_dir_name,
                                           meta_output_dir,
                                            '20_orfs_plots/orfs_gc_contents.txt',
                                           figure_str_to_show_on_html='GC content per genome dispersion')

    html_text += get_html_string_of_result(final_output_dir_name,
                                           meta_output_dir,
                                            '11_final_table/final_orthologs_table.csv',
                                           raw_str_to_show_on_html='Orthologs groups table')

    html_text += get_html_string_of_result(final_output_dir_name,
                                           meta_output_dir,
                                            '11_final_table/phyletic_pattern.fas',
                                           raw_str_to_show_on_html='Phyletic pattern',
                                           additional_text='&nbsp;(Further analyze gain/loss dynamics with <a href="http://gloome.tau.ac.il/" target="_blank">GLOOME</a>)')

    html_text += f'</tbody></table>\n</div>\n'

    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()

    add_closing_html_tags(html_path, run_number)


def edit_failure_html(html_path, run_number, msg):
    if consts.IGNORE_HTML:
        return

    html_text = ''
    try:
        with open(html_path) as f:
            html_text = f.read()
        # The initial file exists (generate by the cgi) so we can read and parse it.
        html_text = html_text.replace('RUNNING', 'FAILED').replace(
            f'{consts.WEBSERVER_NAME} is now processing your request. This page will be automatically updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. ',
            '')
    except FileNotFoundError:
        import logging
        logger = logging.getLogger('main')
        logger.warning(f"Couldn't find html prefix at: {html_path}")

    html_text += f'<div class="container" align="justify" style="{consts.CONTAINER_STYLE}"><h3>\n' \
                 f'<font color="red">{msg}</font><br><br>' \
                 f'Please make sure your input is OK and then try to re-run your job or <a href="mailto:{consts.ADMIN_EMAIL}?subject={consts.WEBSERVER_NAME}%20Run%20Number:%20{run_number}">contact us</a> for further information' \
                 f'</h3></div>\n'

    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()

    add_closing_html_tags(html_path, run_number)


def edit_progress(output_html_path, progress=None, active=True):
    if consts.IGNORE_HTML:
        return

    result = ''
    with open(output_html_path) as f:
        for line in f:
            if 'progress-bar' in line:
                if progress:
                    line = line.split('style')[0]  # <div class="progress-bar ... style="width:0%">\n
                    line += f'style="width:{progress}%">\n'
                if not active:
                    line = line.replace('progress-bar-striped active',
                                        'progress-bar-striped')  # <div class="progress-bar progress-bar-striped active" ...
            result += line

    with open(output_html_path, 'w') as f:
        f.write(result)
