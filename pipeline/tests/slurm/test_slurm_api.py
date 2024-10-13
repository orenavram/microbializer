import os
import sys
import json
import pathlib
import logging
import platform

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from pipeline.flask.slurm_example import submit_job
from pipeline.flask.flask_interface_consts import MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME, NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB, \
    MICROBIALIZER_PROCESSOR_JOB_PREFIX, MICROBIALIZER_PROCESSOR_RESULTS_FILE_NAME, MICROBIALIZER_JOB_TEMPLATE, \
    ARGS_JSON_PATH_KEY, JOB_PARAMETERS_FILE_NAME
from pipeline.flask.SharedConsts import INTERVAL_BETWEEN_LISTENER_SAMPLES


LOG_FILE = '/groups/pupko/yairshimony/test/test_slurm_api/log.txt' if platform.system() == 'Linux' else r'C:\repos\microbializer\pipeline\tests\test_slurm_api_log.txt'
GENOMES_PATH = "/lsweb/pupko/microbializer/user_results/test_slurm_api/4_genomes.zip"


def test_slurm_example(input_path, job_arguments, logger):
    input_path_parent = os.path.dirname(input_path)
    job_unique_id = os.path.basename(input_path_parent)
    results_file_path = f"{input_path_parent}/{MICROBIALIZER_PROCESSOR_RESULTS_FILE_NAME}"
    job_name = f'{MICROBIALIZER_PROCESSOR_JOB_PREFIX}_{job_unique_id}'
    json_parameters_file_path = f"{input_path_parent}/{JOB_PARAMETERS_FILE_NAME}"
    # with open(json_parameters_file_path, 'w') as fp:
    #     json.dump(job_arguments, fp)
    command_to_run = MICROBIALIZER_JOB_TEMPLATE.format(
        sleep_interval=INTERVAL_BETWEEN_LISTENER_SAMPLES,
        args_json_path_key=ARGS_JSON_PATH_KEY,
        args_json_path=json_parameters_file_path,
        results_file_path=results_file_path
    )

    # run the job
    run_parameters = {"queue": MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME,
                      "num_cpus": NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB, "job_name": job_name,
                      "logs_path": input_path_parent, "script_commands": command_to_run, "memory": "10",
                      'logger': logger, 'current_working_directory': input_path_parent}
    logger.debug(f'{run_parameters}')
    return submit_job(**run_parameters)


if __name__ == '__main__':
    job_arguments = {
        "contigs_dir": GENOMES_PATH
    }
    logging.basicConfig(filename=LOG_FILE, level=logging.INFO)
    logger = logging.getLogger(__name__)
    test_slurm_example(GENOMES_PATH, job_arguments, logger)
