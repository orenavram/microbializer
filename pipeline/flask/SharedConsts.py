from pathlib import Path
from enum import Enum


WEBSERVER_DOMAIN = 'microbializer.tau.ac.il'
WEBSERVER_ADDRESS = f'https://{WEBSERVER_DOMAIN}'
OWNER_EMAIL = 'yairshimony@mail.tau.ac.il'

# OUTPUT consts
K_MER_COUNTER_MATRIX_FILE_NAME = Path('CounterMatrixForUI.csv')
RESULTS_FOR_OUTPUT_CLASSIFIED_RAW_FILE_NAME = Path('ResultsForPostProcessClassifiedRaw.csv')
RESULTS_SUMMARY_FILE_NAME = Path('summary_results.txt')
RESULTS_FOR_OUTPUT_UNCLASSIFIED_RAW_FILE_NAME = Path('ResultsForPostProcessUnClassifiedRaw.csv')
TEMP_CLASSIFIED_IDS = Path('TempClassifiedIds.txt')
TEMP_UNCLASSIFIED_IDS = Path('TempUnClassifiedIds.txt')
INPUT_CLASSIFIED_FILE_NAME = Path('classified.fasta')
INPUT_UNCLASSIFIED_FILE_NAME = Path('unclassified.fasta')
INPUT_CLASSIFIED_FILE_NAME_PAIRED = Path('classified_paired.fasta')
INPUT_UNCLASSIFIED_FILE_NAME_PAIRED = Path('unclassified_paired.fasta')
FINAL_OUTPUT_FILE_NAME = Path('FilteredResults.txt.gz')
NEW_CONTAMINATION_FASTA = Path('NewContaminatedSeqs.fasta')
NEW_CONTAMINATION_FASTA_PAIRED = Path('NewContaminatedSeqsPaired.fasta')
FINAL_OUTPUT_FILE_CONTAMINATED_NAME = Path('NewContaminatedSeqs.fasta.gz')
FINAL_OUTPUT_ZIPPED_BOTH_FILES_NEW_CONTAMINATED = Path('NewContaminatedSeqs.tar.gz')
FINAL_OUTPUT_FILE_NAME_PAIRED = Path('FilteredResults_paired.txt.gz')
FINAL_OUTPUT_ZIPPED_BOTH_FILES = Path('FilteredResults_paired_files.tar.gz')
KRAKEN_SUMMARY_RESULTS_FOR_UI_FILE_NAME = Path('summary_stat_UI.json')
GENOME_DOWNLOAD_SUMMARY_RESULTS_FILE_NAME = Path('genome_download_summary.txt')
RANK_KRAKEN_TRANSLATIONS = {'U': 'Unclassified', 'R': 'Root', 'D': 'Domain', 'K': 'Kingdom', 'P': 'Phylum',
                            'C': 'Class', 'O': 'Order', 'F': 'Family', 'G': 'Genus', 'S': 'Species'}

USER_FILE_NAME_ZIP = "genomes.zip"  # names of the files to be downloaded
USER_FILE_NAME_TAR = "genomes.tar.gz"  # names of the files to be downloaded
MAX_NUMBER_PROCESS = 15  # number of process to run in parallel

PATH_TO_OUTPUT_PROCESSOR_SCRIPT = Path("/bioseq/genome_fltr_backend/OutputProcessor.py")
CUSTOM_DB_NAME = 'custom'  # the name of the custom db type

DF_LOADER_CHUCK_SIZE = 1e6
RESULTS_COLUMNS_TO_KEEP = ['is_classified', 'read_name', 'max_specie', 'classified_species', 'read_length',
                           'max_k_mer_p',
                           'all_classified_K_mers', 'split']
SUMMARY_RESULTS_COLUMN_NAMES = ['percentage_of_reads', 'number_of_reads_under', 'number_of_reads_exact', 'rank_code',
                                'ncbi_taxonomyID', 'name']
UNCLASSIFIED_COLUMN_NAME = 'Non Bacterial'
KRAKEN_UNCLASSIFIED_COLUMN_NAME = 'unclassified (taxid 0)'
UNCLASSIFIED_BACTERIA_NAME = 'Unclassified Contamination'
UNCLASSIFIED_BACTERIA_ID = '-1'

# PBS Listener consts
JOB_NUMBER_COL = 'job_id'
JOB_NAME_COL = 'name'
JOB_STATUS_COL = 'current_state'
JOB_ELAPSED_TIME = 'elapsed_time'
JOB_CHANGE_COLS = [JOB_NUMBER_COL, JOB_NAME_COL, JOB_STATUS_COL]
QstatDataColumns = [JOB_NUMBER_COL, 'username', 'queue', JOB_NAME_COL, 'session_id', 'nodes', 'cpus', 'req_mem',
                    'req_time', JOB_STATUS_COL, JOB_ELAPSED_TIME]
SRVER_USERNAME = 'bioseq'
JOB_RUNNING_TIME_LIMIT_IN_HOURS = 10

# Job listener and management function naming
LONG_RUNNING_JOBS_NAME = 'LongRunning'
QUEUE_JOBS_NAME = 'Queue'
NEW_RUNNING_JOBS_NAME = 'NewRunning'
FINISHED_JOBS_NAME = 'Finished'
ERROR_JOBS_NAME = 'Error'
WEIRD_BEHAVIOR_JOB_TO_CHECK = ''
PATH2SAVE_PROCESS_DICT = r'SavedObjects/processes.dict'
PATH2SAVE_PREVIOUS_DF = r'SavedObjects/previous_processes.csv'
PATH2SAVE_WAITING_LIST = r'SavedObjects/waiting.lst'
PATH2SAVE_MONITOR_DATA = r'SavedObjects/monitored_data'

INTERVAL_BETWEEN_CLEANING_THE_PROCESSES_DICT = 24  # in hours
TIME_TO_SAVE_PROCESSES_IN_THE_PROCESSES_DICT = 30  # in days
TIME_TO_KEEP_PROCSES_IDS_FOLDERS = 14  # in days the entire folder of the process

# post processing
POSTPROCESS_JOB_PREFIX = 'PP'
POSTPROCESS_JOB_QUEUE_NAME = 'lifesciweb'
NUBMER_OF_CPUS_POSTPROCESS_JOB = '1'
NEW_CONTAMINATION_PAIRED_CREATION_LINE = """
cp "$original_classified_data_paired" "$output_path_new_classified_data_temp"\n
cat "$output_path_new_classified_data_temp" | seqkit grep -f "$output_pathTemp" -v -o "$output_path_new_classified_data_paired"
"""
FILTER_ORIGINAL_CLASSIFIED_RESULTS = 'cat "$original_classified_data_paired" | seqkit grep -f "$output_pathTemp"  -o "$Temp_new_unclassified_seqs_paired"\n'
COMBINE_RESULTS_PAIRED = 'cat "$Temp_new_unclassified_seqs_paired" "$original_unclassified_data_paired" > "$output_pathPaired"\n'
GZIP_LINE_PAIRED = '''
cd "{path_to_process_folder}"\n
tar cvzf "{tar_result_file_name}" "{output_one_file_name}" "{output_paired_file_name}"\n
'''
GZIP_LINE_PAIRED_NEW_CONTAMINATED = '''
cd "{path_to_process_folder}"\n
tar cvzf "{tar_result_file_name}" "{output_one_file_name}" "{output_paired_file_name}"\n
'''
GZIP_LINE_ONE_FILE = 'gzip -f "$output_path"\n'
GZIP_LINE_ONE_FILE_CONTAMINATED = 'gzip -f "$output_path_new_classified_data"\n'
POST_PROCESS_COMMAND_TEMPLATE = '''
#!/bin/bash          
#PBS -S /bin/bash
#PBS -r y
#PBS -q {queue_name}
#PBS -l ncpus={cpu_number}
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N {job_name}
#PBS -e {error_files_path}
#PBS -o {output_files_path}

source /powerapps/share/miniconda3-4.7.12/etc/profile.d/conda.sh
conda activate NGScleaner

sleep {sleep_interval}

original_unclassified_data="{path_to_original_unclassified_data}"
original_classified_data="{path_to_original_classified_data}"
original_unclassified_data_paired="{path_to_original_unclassified_data_paired}"
original_classified_data_paired="{path_to_original_classified_data_paired}"
input_path="{path_to_classified_results}"
output_path="{path_to_final_result_file}"
output_pathTemp="{path_to_temp_file}"
output_pathPaired="{path_to_final_result_file_paired}"
Temp_new_unclassified_seqs="{path_to_temp_unclassified_file}"
Temp_new_unclassified_seqs_paired="{path_to_temp_unclassified_file_paired}"
string='{species_to_filter_on}'
output_path_new_classified_data='{path_to_new_contamination_seqs}'
output_path_new_classified_data_temp='{path_to_temp_new_contamination_seqs}'
output_path_new_classified_data_paired='{path_to_new_contamination_seqs_paired}'

# filter kraken results by query name and threshold
cat "$input_path" | awk -F "," '{{split(var,parts,","); for (i in parts) dict[parts[i]]; if ($6 <= {classification_threshold} || !($3 in dict)) print }}' var="${{string}}" | awk -F "," 'NR!=1 {{print $2}}' > "$output_pathTemp"

# filter original classified results
cat "$original_classified_data" | seqkit grep -f "$output_pathTemp"  -o "$Temp_new_unclassified_seqs"
{filter_original_classified_results_line}

# create a copy of original classified data without the seqs the user defined as NOT contamination
cp "$original_classified_data" "$output_path_new_classified_data_temp"
cat "$output_path_new_classified_data_temp" | seqkit grep -f "$output_pathTemp" -v -o "$output_path_new_classified_data"
{create_new_contamination_paired_line}

# combine original unfiltered input with newly unclassified results
cat "$Temp_new_unclassified_seqs" "$original_unclassified_data" > "$output_path"
{combine_results_paired_line}

sleep 30

{gzip_line} 
{gzip_line_new_contaminated}

#rm "$output_pathTemp"
#rm "$Temp_new_unclassified_seqs"
'''

SEPERATOR_FOR_MONITOR_DF = '###'


class State(Enum):
    Running = 1
    Finished = 2
    Crashed = 3
    Waiting = 4
    Init = 5
    Queue = 6


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
        elif state.value == State.Init.value:
            if job_name:
                return f'Microbializer {job_name} - Job Initialized'
            return f'Microbializer - Job Initialized'
        elif state.value == State.Running.value:
            if job_name:
                return f'Microbializer {job_name} - Job Running'
            return f'Microbializer - Job Running'
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
    static_folder_path = 'gifs/'
    states_gifs_dict = {
        State.Running: {
            "background": "#ffffff",
            "gif_id": "D5GyCFkInbJlu"
        },
        State.Finished: {
            "background": "#1674d2",
            "gif_id": "TvLuZ00OIADoQ"
        },
        State.Crashed: {
            "background": "#1674d2",
            "gif_id": "TvLuZ00OIADoQ"
        },
        State.Waiting: {
            "background": "#1674d2",
            "gif_id": "TvLuZ00OIADoQ"
        },
        State.Init: {
            "background": "#1674d2",
            "gif_id": "TvLuZ00OIADoQ"
        },
        State.Queue: {
            "background": "#1674d2",
            "gif_id": "TvLuZ00OIADoQ"
        }
    }

    states_text_dict = {
        State.Running: "Your process is running...",
        State.Finished: "Your process finished... Redirecting to results page.",  # TODO is needed??
        State.Crashed: "Your process failed.\nWe suggest you rerun the process.",  # TODO finish
        State.Waiting: "We are currently running other processes.\nYour process will start soon.",
        State.Init: "The process is initialized and will shortly enter the next stage.",
        State.Queue: "Job is queued",
    }

    global allowed_files_str  # todo: Edo, do we have to use a global var?
    ALLOWED_EXTENSIONS = {'zip', 'gz'}
    allowed_files_str = ', '.join(ALLOWED_EXTENSIONS)  # better to path string than list

    class UI_Errors(Enum):
        UNKNOWN_PROCESS_ID = 'The provided process id does not exist'
        INVALID_EXPORT_PARAMS = 'invalid paramters for export'
        POSTPROCESS_CRASH = 'can\'t postprocess'
        INVALID_MAIL = 'invalid mail'
        CANT_ADD_PROCESS = 'can\'t add search process'
        INVALID_FILE_EXTENTION = f'invalid file extenstion, please upload a file of: {allowed_files_str}'
        CORRUPTED_FILE = f'please upload a valid fasta file, the file you uploaded is corrupted'
        INVALID_FILES_NUMBER = f'please insert one or two files'
        EXPORT_FILE_UNAVAILABLE = f'failed to export file, try to rerun the file'
        HISTOGRAM_DATA_IS_NULL = f'cannot find the data for histogram, make sure the process id is correct and finished'
        ORTHOLOGOUS_DATA_IS_NULL = f'cannot find the data for orthologous, make sure the process id is correct and finished'
        NEWICK_DATA_IS_NULL = f'cannot find the data for newick tree, make sure the process id is correct and finished'
        EXPORT_FILE_CONTAMINATED_UNAVAILABLE = f'failed to export the contaminated file, this may be because the process run before the update of the contaminated file'
        PAGE_NOT_FOUND = 'The requested page does not exist'
        INVALID_SPECIES_LIST = 'Some of the species inserted to the custom DB are invalid'
        JOB_CRASHED = 'Your processed crashed... Make sure your input is valid'
        RESULTS_DF_IS_NONE = 'The results cant be found, can be caused by several things, please contact us'
        UNKNOWN_ACTION = 'Please use send a valid action'
        NO_ACTION = 'No action was added, please choose'
        ALL_FILES_NOT_CREATE = 'Cannot find the required zip files for all results'
        FILE_NOT_FOUND = 'Cannot find the required file'

    ERROR_CONTACT_INFO = f'For more information, or any other inquiries, please contact {OWNER_EMAIL}'

    PROCESS_INFO_KR = f"We are processing your request. This may take an hour for a small number of genomes and up to several days for a large one. This link is valid for {TIME_TO_KEEP_PROCSES_IDS_FOLDERS} days. if an email address was provided, a link will be sent upon analysis completion."

    TEXT_TO_RELOAD_HTML = "update"  # empty string not allowed! will cause massive bug!
    FETCH_UPDATE_INTERVAL_HTML_SEC = 15  # should be greater than 5 (keep-alive timeout mod_wsgi).

    #  About page text
    HELP_TEXT_ABOUT_LIST = [
        """Here we present the GenomeFLTR web server. The web server was developed to easily filter genomic reads at a click of a mouse. The server automatically updates the databases and provides a simple and interactive GUI (graphical user interface) to personalize the user output. Visual and textual results that are ready for publication or further analysis are provided."""]
    HELP_TEXT_JOB_NAME = """The job name as inserted by the user"""
    HELP_TEXT_TAXA_DOWNLOAD = """The referenced genome downloaded to compare against"""
