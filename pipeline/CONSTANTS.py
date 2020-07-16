#!/powerapps/share/centos7/python-anaconda3.6.5/bin/python

#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

OWNER_EMAIL = 'orenavram@gmail.com'

# general vars
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'
RELOAD_INTERVAL = 5
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

# relevant modules
GCC = 'gcc/gcc-7.3.0'
MCL = 'MCL-edge/mcl-14-137'
MAFFT = 'mafft/mafft-7.407'
RAXML = 'raXML'
PRODIGAL = 'prodigal/prodigal-2.6.3'
MMSEQS = 'MMseqs2/June2020'

WEBSERVER_NAME = 'M1CR0B1AL1Z3R'
WEBSERVER_URL = 'https://microbializer.tau.ac.il'
WEBSERVER_TITLE = 'A web server for analyzing bacterial genomics data. Easily.'

WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'microbializer')
WEBSERVER_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'microbializer')
WEBSERVER_HTML_DIR = '/data/www/html/microbializer'

WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')

MAIN_SCRIPT = '/bioseq/microbializer/pipeline/microbializer_pipeline.py'
SUBMISSIONS_LOG = '/bioseq/microbializer/submissions_log.txt'
EMAIL_FILE_NAME = 'email.txt'

#path to example data
EXAMPLE_DATA = os.path.join(WEBSERVER_HTML_DIR, 'example_data.zip')

WEBSERVER_JUMBOTRON = f'&nbsp;&nbsp;&nbsp;&nbsp;<span id="server-title">{WEBSERVER_NAME}</span>&nbsp;&nbsp;&nbsp;&nbsp;<span id="sub-title">{WEBSERVER_TITLE}</span>'

PROGRESS_BAR_TAG = '''<div class="progress">
        <div class="progress-bar progress-bar-striped active" role="progressbar" aria-valuemin="0" aria-valuemax="100" style="width:5%">
ORFs detection ; Searching homologs ; Clustering homologs ; Aligning homologs ; Core proteome inferrence ; Tree reconstruction ; Additional statistics
        </div>
    </div>'''

CONTAINER_WIDTH = 'width: 850px'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}'
