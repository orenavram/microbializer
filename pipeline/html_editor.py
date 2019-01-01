import os


def add_closing_html_tags(html_path, CONSTS, run_number):
    with open(html_path, 'a') as f:
        f.write(
            f'<hr>\n<h4 class=footer><p align=\'center\'>Questions and comments are welcome! Please ' \
            f'<span class="admin_link">' \
            f'<a href="mailto:{CONSTS.ADMIN_EMAIL}?subject=ASAP%20Run%20Number%20{run_number}">contact us</a>' \
            f'</span></p></h4>\n' \
            f'<div id="bottom_links" align="center"><span class="bottom_link">' \
            f'<a href="{CONSTS.WEBSERVER_URL}/" target="_blank">Home</a>' \
            f'&nbsp;|&nbsp<a href="{CONSTS.WEBSERVER_URL}/overview.html" target="_blank">Overview</a>\n' \
            f'</span>\n' \
            f'<br><br><br>\n</body>\n</html>\n')
        f.flush()

    # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
    from time import sleep
    sleep(2 * CONSTS.RELOAD_INTERVAL)
    with open(html_path) as f:
        html_content = f.read()
    html_content = html_content.replace(CONSTS.RELOAD_TAGS, '')
    with open(html_path, 'w') as f:
        f.write(html_content)


def get_html_string_of_restult(final_output_dir_name, meta_output_dir, end_of_str, figure_str_to_show_on_html='', raw_str_to_show_on_html='raw data'):
    result = '<tr><td>'
    raw_file_suffix = os.path.join(final_output_dir_name, end_of_str)
    if os.path.exists(os.path.join(meta_output_dir, raw_file_suffix)):
        if figure_str_to_show_on_html:
            result += f'<a href="{raw_file_suffix.replace("txt", "png")}" target="_blank">{figure_str_to_show_on_html}</a> ; ('
        result += f'<a href="{raw_file_suffix}" target="_blank">{raw_str_to_show_on_html}</a>'
        if figure_str_to_show_on_html:
            result += ')\n'
        result += f'<br></td></tr>'

    return result


def edit_success_html(html_path, meta_output_dir, final_output_dir_name, run_number, CONSTS):
    html_text = ''
    try:
        with open(html_path) as f:
            html_text = f.read()
        # The initial file exists (generate by the cgi) so we can read and parse it.
        html_text = html_text.replace('RUNNING', 'FINISHED').replace(f'{CONSTS.WEBSERVER_NAME} is now processing your request. This page will be automatically updated every {CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. ', '')
    except FileNotFoundError:
        import logging
        logger = logging.getLogger('main')
        logger.warning(f"Couldn't find html prefix at: {html_path}")

    html_text += f'<div class="container" style="{CONSTS.CONTAINER_STYLE}">\n' \
        f'<h2>RESULTS:<h2>'\
        f'<h3><b><a href=\'{CONSTS.WEBSERVER_NAME}_outputs.zip\' target=\'_blank\'>Download zipped full results (textual & visual)</a></b></h3>' \
        f'<table class="table">\n' \
        f'<thead>\n' \
        f'<tr><th><h3>Quick access to selected results:</h3></th></tr>\n' \
        f'</thead>\n' \
        f'<tbody>'

    html_text += get_html_string_of_restult(final_output_dir_name,
                                            meta_output_dir,
                                            '11_final_table/final_orthologs_table.csv',
                                            raw_str_to_show_on_html='Orthologs groups table')

    html_text += get_html_string_of_restult(final_output_dir_name,
                                            meta_output_dir,
                                            '13_groups_sizes_frequency/groups_sizes_frequency.txt',
                                            figure_str_to_show_on_html='Orthologs groups size dispersion')

    html_text += get_html_string_of_restult(final_output_dir_name,
                                            meta_output_dir,
                                            '14_orfs_statistics/orfs_counts.txt',
                                            figure_str_to_show_on_html='ORFs per genome dispersion')

    html_text += get_html_string_of_restult(final_output_dir_name,
                                            meta_output_dir,
                                            '14_orfs_statistics/orfs_gc_contents.txt',
                                            figure_str_to_show_on_html='GC content per genome dispersion')

    html_text += get_html_string_of_restult(final_output_dir_name,
                                            meta_output_dir,
                                            '17_species_phylogeny/species_tree.txt',
                                            figure_str_to_show_on_html='Species tree')

    html_text += f'</tbody></table>\n' \
        f'</div>\n'

    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()

    add_closing_html_tags(html_path, CONSTS, run_number)



def edit_failure_html(html_path, run_number, msg, CONSTS):
    html_text = ''
    try:
        with open(html_path) as f:
            html_text = f.read()
        # The initial file exists (generate by the cgi) so we can read and parse it.
        html_text = html_text.replace('RUNNING', 'FAILED').replace(f'{CONSTS.WEBSERVER_NAME} is now processing your request. This page will be automatically updated every {CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. ', '')
    except FileNotFoundError:
        import logging
        logger = logging.getLogger('main')
        logger.warning(f"Couldn't find html prefix at: {html_path}")

    html_text +=f'<br><br><br>\n' \
                f'<center><h2>\n' \
                f'<font color="red">{msg}</font><br><br>' \
                f'Please try to re-run your job or <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject=ASAP%20Run%20Number%20{run_number}">contact us</a> for further information' \
                f'</h2></center>\n' \
                f'<br><br>\n' \
                f'</body>\n</html>\n'

    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()

    add_closing_html_tags(html_path, CONSTS, run_number)

