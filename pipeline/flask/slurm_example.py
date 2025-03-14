import requests
import os
import api_secrets
import flask_interface_consts

# Base URL for authentication and token generation
base_url_auth = 'https://slurmtron.tau.ac.il'
generate_token_url = f"{base_url_auth}/slurmapi/generate-token/"
# Base URL for job submission
base_url = f"{base_url_auth}/slurmrestd"
# Job submission URL
job_url = f'{base_url}/slurm/v0.0.40/job/submit'

# User credentials
current_user = "microbi"


def get_api_token(username, api_key):
    """
    Retrieves a JWT token for SLURM REST API access for powerslurm cluster.

    Parameters:
    username (str): The username of the user requesting the token.
    api_key (str): The API key provided by the HPC team.

    Returns:
    str: The API token if the request is successful.

    Raises:
    Exception: If the request fails with a non-200 status code.
    """

    generate_token_url = 'https://slurmtron.tau.ac.il/slurmapi/generate-token/'

    payload = {
        'username': username,
        'api_key': api_key
    }

    response = requests.post(generate_token_url, data=payload)

    if response.status_code == 200:
        # Extracting the token from the JSON response
        return response.json()['SlurmJWT']
    else:
        raise Exception(f"Error: {response.status_code}, {response.text}")


def submit_job(script_commands, job_name, logs_path, logger, current_working_directory="/tmp/"):
    logger.info(f'in submit_job, for {job_name}')
    # Authorization headers with the obtained token
    headers = {
        'X-SLURM-USER-NAME': current_user,
        'X-SLURM-USER-TOKEN': get_api_token(current_user, api_secrets.API_KEY)
    }

    slurm_output_file = os.path.join(logs_path, 'main_%j.out')
    slurm_error_file = os.path.join(logs_path, 'main_%j.err')

    slurm_script_path = os.path.join(logs_path, 'job.slurm')
    with open(slurm_script_path, 'w') as f:
        slurm_header = flask_interface_consts.MICROBIALIZER_JOB_HEADER_TEMPLATE.format(
            output_file=slurm_output_file, error_file=slurm_error_file)
        shebang, commands = script_commands.split('\n', 1)
        full_slurm_script = f"{shebang}\n{slurm_header}\n{commands}"
        f.write(full_slurm_script)
        logger.info(f"Wrote slurm script to {slurm_script_path}")

    # Job submission request
    jobs_request = requests.post(
        job_url,
        headers=headers,
        json={
            # Example job script
            "script": script_commands,
            "job": {
                "partition": flask_interface_consts.MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME,
                "tasks": 1,
                "name": job_name,
                # "account": "< account_name >",
                "nodes": "1",
                # "allocation_node_list": "compute-0-299",
                "cpus_per_task": int(flask_interface_consts.NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB),
                "memory_per_node": {
                    "number": str(flask_interface_consts.MICROBIALIZER_MAIN_JOB_MEMORY),
                    "set": True,
                    "infinite": False
                },
                # we should pass the time_limit in minutes
                "time_limit": flask_interface_consts.MICROBIALIZER_MAIN_JOB_TIME_LIMIT_IN_HOURS * 60,
                # Full path to your error/output file.
                "standard_output": slurm_output_file,
                "standard_error": slurm_error_file,
                "current_working_directory": current_working_directory,
                # Environment modules (module load) should not be used directly under the script parameter. Instead, set all necessary environment variables under the environment parameter.
                "environment": [
                    "PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/root/bin:",
                    # "PATH=/powerapps/share/rocky8/gcc-13.1.0/bin:/lsweb/pupko/microbializer/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/root/bin:", #/a/home/cc/lifesci/microbializer/.pyenv/shims:/a/home/cc/lifesci/microbializer/.pyenv/bin
                    # "LD_LIBRARY_PATH=/lib/:/lib64/:/usr/local/lib:/powerapps/share/rocky8/gcc-13.1.0/lib64:/powerapps/share/rocky8/gcc-13.1.0/lib:/powerapps/src/rocky8/gcc-13.1.0/isl/lib:/powerapps/src/rocky8/gcc-13.1.0/mpfr-4.2.0/lib:/powerapps/src/rocky8/gcc-13.1.0/gmp-6.2.1/lib"
                ],
            },
        }
    )

    # Processing the job submission result
    jobs_result = jobs_request.json()
    if 'result' in jobs_result:
        jobs_request = jobs_result['result']
        logger.info(jobs_result)
        for key, value in jobs_result.items():
            if 'job_id' in key:
                logger.info(f'submitted job_id {value}')
                return value
    else:
        logger.error(f'failed to submit job {job_name} return from job = {jobs_result}')
