from enum import Enum


WEBSERVER_DOMAIN = 'microbializer.tau.ac.il'
WEBSERVER_ADDRESS = f'https://{WEBSERVER_DOMAIN}'
OWNER_EMAIL = 'yairshimony@mail.tau.ac.il'

ALLOWED_EXTENSIONS = ['zip', 'gz']
USER_FILE_NAME_ZIP = "genomes.zip"  # names of the files to be downloaded
USER_FILE_NAME_TAR = "genomes.tar.gz"  # names of the files to be downloaded

TIME_TO_KEEP_PROCSES_IDS_FOLDERS = 14  # in days the entire folder of the process


class State(Enum):
    Running = 1
    Finished = 2
    Crashed = 3
    Init = 5


class EMAIL_CONSTS:
    def create_title(state: State, job_name):
        if state.value == State.Finished.value:
            if job_name:
                return f'Microbializer {job_name} - Job Finished'
            return f'Microbializer - Job Finished'
        elif state.value == State.Crashed.value:
            if job_name:
                return f'Microbializer {job_name} - Job Crashed'
            return f'Microbializer - Job Crashed'
        else:
            return f'unknown state in create_title at EMAIL_CONSTS'

    CONTENT_PROCESS_CRASHED = f"Thank you for using Microbializer.\n" \
                              f"We are sorry for the inconvenience, but the process crashed.\n" \
                              f"Please look at: {WEBSERVER_ADDRESS}/error_from_job/{{process_id}} for the error details and verify your input.\n" \
                              f"For more help contact: {OWNER_EMAIL}"

    CONTENT_PROCESS_FINISHED = f"Thank you for using Microbializer.\n" \
                               f"Your results visual summary is at: {WEBSERVER_ADDRESS}/results/{{process_id}}\n" \
                               f"Your downloadable results are at: {WEBSERVER_ADDRESS}/download_page/{{process_id}}\n" \
                               f"Please note that the results are available for {str(TIME_TO_KEEP_PROCSES_IDS_FOLDERS)} days, " \
                               f"and will be deleted after, so remember to download them (we cannot keep them for a longer time due to storage limitations).\n\n" \
                               f"Please remember to cite us in your work (citation info is at: {WEBSERVER_ADDRESS}/about).\n" \
                               f"For more help contact: {OWNER_EMAIL}"

    SUBMITTED_TITLE = '''Microbializer {job_name} - Job Submitted'''
    SUBMITTED_CONTENT = f"Thank you for using Microbializer.\n" \
                        f"Your job has been submitted, you can check its status at: {WEBSERVER_ADDRESS}/process_state/{{process_id}}\n" \
                        f"An update will be sent upon completion."


class UI_CONSTS:
    states_gifs_dict = {
        State.Running: {
            "background": "#ffffff",
            "gif_id": "D5GyCFkInbJlu"
        },
        State.Init: {
            "background": "#1674d2",
            "gif_id": "TvLuZ00OIADoQ"
        }
    }

    states_text_dict = {
        State.Running: "Your process is running...",
        State.Init: "The job is submitted and will soon enter the next stage.\nIt might take a while, "
                    "depending on the number of other jobs currently running.",
    }

    class UI_Errors(Enum):
        UNKNOWN_PROCESS_ID = 'The provided process id does not exist'
        INVALID_MAIL = 'invalid mail'
        INVALID_FILE_EXTENTION = f'invalid file extenstion, please upload a file of: {", ".join(ALLOWED_EXTENSIONS)}'
        CORRUPTED_FILE = f'please upload a valid fasta file, the file you uploaded is corrupted'
        INVALID_FILES_NUMBER = f'please insert one or two files'
        HISTOGRAM_DATA_IS_NULL = f'cannot find the data for histogram, make sure the process id is correct and finished'
        ORTHOLOGOUS_DATA_IS_NULL = f'cannot find the data for orthologous, make sure the process id is correct and finished'
        NEWICK_DATA_IS_NULL = f'cannot find the data for newick tree, make sure the process id is correct and finished'
        PAGE_NOT_FOUND = 'The requested page does not exist'
        JOB_CRASHED = 'Your processed crashed... Make sure your input is valid'
        FILE_NOT_FOUND = 'Cannot find the required file'

    ERROR_CONTACT_INFO = f'For more information, or any other inquiries, please contact {OWNER_EMAIL}'

    PROCESS_INFO_KR = f"We are processing your request. This may take an hour for a small number of genomes and up to several days for a large one. This link is valid for {TIME_TO_KEEP_PROCSES_IDS_FOLDERS} days. if an email address was provided, a link will be sent upon analysis completion."

    TEXT_TO_RELOAD_HTML = "update"  # empty string not allowed! will cause massive bug!
    FETCH_UPDATE_INTERVAL_HTML_SEC = 15  # should be greater than 5 (keep-alive timeout mod_wsgi).
