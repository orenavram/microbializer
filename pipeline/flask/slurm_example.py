import requests
import os
import api_secrets
import flask_interface_consts

# Job submission endpoint
job_submit_url = "https://saw.tau.ac.il/slurmapi/job/submit/"

# User credentials
current_user = "microbi"


def submit_job(script_commands, job_name, logs_path, logger):
    logger.info(f'in submit_job, for {job_name}')

    # Authorization headers
    headers = {
        "Authorization": f"Bearer {api_secrets.API_KEY}",
        "X-USERNAME": current_user,
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
    payload = {
        "script": script_commands,
        "partition": flask_interface_consts.MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME,
        "tasks": 1,
        "name": job_name,
        "account": flask_interface_consts.MICROBIALIZER_PROCESSOR_JOB_ACCOUNT_NAME,
        "nodes": 1,
        "qos": flask_interface_consts.MICROBIALIZER_PROCESSOR_JOB_QOS,
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
        "current_working_directory": str(logs_path),
        # Environment modules (module load) should not be used directly under the script parameter. Instead, set all necessary environment variables under the environment parameter.
        "environment": [
            "PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/root/bin:",
            # "PATH=/powerapps/share/rocky8/gcc-13.1.0/bin:/lsweb/pupko/microbializer/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/root/bin:", #/a/home/cc/lifesci/microbializer/.pyenv/shims:/a/home/cc/lifesci/microbializer/.pyenv/bin
            # "LD_LIBRARY_PATH=/lib/:/lib64/:/usr/local/lib:/powerapps/share/rocky8/gcc-13.1.0/lib64:/powerapps/share/rocky8/gcc-13.1.0/lib:/powerapps/src/rocky8/gcc-13.1.0/isl/lib:/powerapps/src/rocky8/gcc-13.1.0/mpfr-4.2.0/lib:/powerapps/src/rocky8/gcc-13.1.0/gmp-6.2.1/lib"
        ],
    }

    jobs_request = requests.post(job_submit_url, headers=headers, json=payload)
    jobs_request.raise_for_status()
    jobs_result = jobs_request.json()

    if 'job_id' in jobs_result:
        job_id = str(jobs_result['job_id'])
        logger.info(f'submitted job_id {job_id}')
        return job_id
    else:
        logger.error(f'failed to submit job {job_name}. return from job submission = {jobs_result}')
