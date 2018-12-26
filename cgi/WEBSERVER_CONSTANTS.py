#!/data/shared/python/anaconda3-5.1.0/bin/python3.6

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

# general paths
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'

# external programs
# IMPORTANT: one must run the command: setenv PATH "/bioseq/Programs/MAFFT_7.222/installation/bin:${PATH}" ahead of this mafft command so all components will be found...
MAFFT_v7_222 = '/bioseq/Programs/MAFFT_7.222/installation/bin/mafft' # v7.222

WEBSERVER_NAME = 'M1CROB1AL1Z3R'
WEBSERVER_URL = 'https://microbializer.tau.ac.il'
#WEBSERVER_LOG = '/bioseq/microbializer/WEBSERVER_runs.log'

RELOAD_INTERVAL = 30
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

OWNER_EMAIL = 'orenavram@gmail.com'

WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'microbializer')
WEBSERVER_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'microbializer')
WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')
WEBSERVER_HTML_DIR = '/data/www/html/microbializer'
WEBSERVER_EXEC = '/bioseq/microbializer/pipeline'

MAIN_SCRIPT = os.path.join(WEBSERVER_EXEC, 'microbializer_pipeline.py')

#path to example data
EXAMPLE_DATA = os.path.join(WEBSERVER_HTML_DIR, 'example_data.tar.gz')
# EXAMPLE_FILE_RUN1_R2 = os.path.join(WEBSERVER_HTML_DIR, 'example', 'run1', '242_R2.fastq')
# EXAMPLE_FILE_RUN2_R1 = os.path.join(WEBSERVER_HTML_DIR, 'example', 'run2', '242_R1.fastq')
# EXAMPLE_FILE_RUN2_R2 = os.path.join(WEBSERVER_HTML_DIR, 'example', 'run2', '242_R2.fastq')

WEBSERVER_TITLE = 'A web server for analyzing bacterial genomics data. Easily.'

#path to IMGT reference library
WEBSERVER_JUMBOTRON = f'&nbsp;&nbsp;&nbsp;&nbsp;<span id="server-title">{WEBSERVER_NAME}</span>&nbsp;&nbsp;&nbsp;&nbsp;<span id="sub-title">{WEBSERVER_TITLE}</span>'