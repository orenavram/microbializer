#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os.path

TEST = True
USE_CONDA = True
IGNORE_HTML = True
CLEAN_OUTPUTS_AFTER_RUN = False
LOG_IN_SEPARATE_FILES = True

# mmseqs command work only on specific machines and this queue navigates only to them
MMSEQS_REQUIRED_MEMORY = '60gb'
QUEUE_FOR_MMSEQS_COMMANDS = 'pupkolab'
USE_DIFFERENT_QUEUE_FOR_MMSEQS = True

# logging consts
JOB_NAME_ENVIRONMENT_VARIABLE = 'PBS_JOBNAME'
JOB_ID_ENVIRONMENT_VARIABLE = 'PBS_JOBID'
LOG_MESSAGE_FORMAT = '%(asctime)s %(levelname)s %(filename)s:%(lineno)d %(message)s'

PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer' if TEST else '/bioseq/microbializer'
OWNER_EMAIL = 'yairshimony@mail.tau.ac.il' if TEST else 'orenavram@gmail.com'

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

# general vars
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'
RELOAD_INTERVAL = 5
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'
CSV_DELIMITER = ','

# relevant modules
GCC = 'gcc/gcc-6.2.0'
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
