import requests
import api_secrets

current_user = "microbi"
base_url = "https://slurmtron.tau.ac.il" # slurprod
generate_token_url = f"{base_url}/slurmapi/generate-token/"
slurmrestd_url = f"{base_url}/slurmrestd"
ACCOUNT_NAME = "pupkoweb-users"

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


def get_jobs(account=None, cluster=None, logger=None):
    #print(f'in get_jobs()')
    #logger.info(f'in get_jobs()')
    params = {}
    if account:
        params['account'] = account
    if cluster:
        params['cluster'] = cluster
    
    headers = {
        'X-SLURM-USER-NAME': current_user,
        'X-SLURM-USER-TOKEN': get_api_token(current_user, api_secrets.API_KEY)
    }
    
    # Sending the GET request to the Slurm REST API
    response = requests.get(f"{slurmrestd_url}/slurmdb/v0.0.40/jobs", headers=headers, params=params)
    # Check if the request was successful
    if response.status_code == 200:
        data = response.json()
        jobs = data.get('jobs', [])
        #logger.info(f"Found {len(jobs)} jobs.")
        #for job in jobs:
        #    logger.info(f"Job ID: {job['job_id']}, Job Name: {job['name']} account: {job['account']}")
    else:
        logger.error(f"Failed to retrieve jobs. Status Code: {response.status_code}, Message: {response.text}")
    return jobs
 
