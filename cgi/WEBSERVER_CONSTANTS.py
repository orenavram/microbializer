#!/data/shared/python/anaconda3-5.1.0/bin/python3.6

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

OWNER_EMAIL = 'orenavram@gmail.com'

# general vars
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'
RELOAD_INTERVAL = 30
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

# external programs
# IMPORTANT: one must run the command: setenv PATH "/bioseq/Programs/MAFFT_7.222/installation/bin:${PATH}" ahead of this mafft command so all components will be found...
MAFFT_v7_222 = '/bioseq/Programs/MAFFT_7.222/installation/bin/mafft' # v7.222


WEBSERVER_NAME = 'M1CR0B1AL1Z3R'
WEBSERVER_URL = 'https://microbializer.tau.ac.il'
WEBSERVER_TITLE = 'A webserver for analyzing bacterial genomics data. Easily.'

WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'microbializer')
WEBSERVER_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'microbializer')
WEBSERVER_HTML_DIR = '/data/www/html/microbializer'

WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')

MAIN_SCRIPT = os.path.join('/bioseq/microbializer/pipeline/microbializer_pipeline.py')

#path to example data
EXAMPLE_DATA = os.path.join(WEBSERVER_HTML_DIR, 'example_data.tar.gz')

WEBSERVER_JUMBOTRON = f'&nbsp;&nbsp;&nbsp;&nbsp;<span id="server-title">{WEBSERVER_NAME}</span>&nbsp;&nbsp;&nbsp;&nbsp;<span id="sub-title">{WEBSERVER_TITLE}</span>'

CONTAINER_WIDTH = 'width: 850px'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}'