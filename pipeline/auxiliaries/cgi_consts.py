# Consts used by CGI and old frontend

import os
from .flask_interface_consts import WEBSERVER_NAME
from .consts import PROJECT_ROOT_DIR


SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'
RELOAD_INTERVAL = 5
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

WEBSERVER_URL = 'https://microbializer.tau.ac.il'
WEBSERVER_TITLE = 'A web server for analyzing bacterial genomics data. Easily.'

WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'microbializer')
WEBSERVER_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'microbializer')
WEBSERVER_HTML_DIR = '/data/www/html/microbializer'

WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')

Q_SUBMITTER_PATH = os.path.join(PROJECT_ROOT_DIR, 'pipeline/auxiliaries/q_submitter_power.py')
MAIN_SCRIPT = '/bioseq/microbializer/pipeline/main.py'
SUBMISSIONS_LOG = '/bioseq/microbializer/submissions_log.txt'
EMAIL_FILE_NAME = 'email.txt'
CGI_DEBUG_FILE_NAME = 'cgi_debug.txt'
RESULT_WEBPAGE_NAME = 'result.html'
EXAMPLE_DATA_FILE_NAME = 'example_data.zip'

# path to example data
EXAMPLE_DATA = os.path.join(WEBSERVER_HTML_DIR, EXAMPLE_DATA_FILE_NAME)

WEBSERVER_JUMBOTRON = f'&nbsp;&nbsp;&nbsp;&nbsp;<span id="server-title">{WEBSERVER_NAME}</span>&nbsp;&nbsp;&nbsp;&nbsp;<span id="sub-title">{WEBSERVER_TITLE}</span>'

PROGRESS_BAR_TAG = '''<div class="progress">
        <div class="progress-bar progress-bar-striped" role="progressbar" aria-valuemin="0" aria-valuemax="100" style="width:0%">
ORFs detection ; Searching homologs ; Clustering homologs ; Aligning homologs ; Core proteome inferrence ; Tree reconstruction ; Additional statistics
        </div>
    </div>'''

CONTAINER_WIDTH = 'width: 850px'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}'