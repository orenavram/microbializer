def add_closing_html_tags(html_path, CONSTS, run_number):
    with open(html_path, 'a') as f:
        f.write(
            f'<br><br><br>\n<hr>\n<h4 class=footer><p align=\'center\'>Questions and comments are welcome! Please ' \
            f'<span class="admin_link">' \
            f'<a href="mailto:{CONSTS.ADMIN_EMAIL}?subject=ASAP%20Run%20Number%20{run_number}">contact us</a>' \
            f'</span></p></h4>\n' \
            f'<div id="bottom_links" align="center"><span class="bottom_link">' \
            f'<a href="{CONSTS.WEBSERVER_URL}/" target="_blank">Home</a>' \
            f'&nbsp;|&nbsp<a href="{CONSTS.WEBSERVER_URL}/overview.html" target="_blank">Overview</a>\n' \
            f'</span>\n' \
            f'<br><br><br>\n</body>\n</html>\n')
        f.flush()


def edit_success_html(html_path, run_number, remote_run, CONSTS):
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

    html_text += f'<div class="container" align="center" style="{CONSTS.CONTAINER_STYLE}">\n' \
        f'<br><br><center><h2>RESULTS:<h2>'\
        f'<a href=\'{CONSTS.WEBSERVER_NAME}_outputs.zip\' target=\'_blank\'><h3><b>Download zipped full results</b></h3></a>' \
        f'</center><br>\n' \
        f'<table class="table">\n' \
        f'<thead>\n' \
        f'<tr><th class="text-center">Analysis Plots</th></tr>\n' \
        f'</thead>\n' \
        f'<tr><td><a href=""></a></td></tr>\n' \
        f'</table>\n' \
        f'</div>\n'

    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()

    add_closing_html_tags(html_path, CONSTS, run_number)



def edit_failure_html(html_path, run_number, msg, remote_run, CONSTS):
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

