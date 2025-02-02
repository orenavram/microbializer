import os
import pandas as pd
import logging


SCRIPT_DIR = Path(__file__).resolve().parent
MAX_PARALLEL_JOBS = 50


def main(logger, max_parallel_jobs, putative_orthologs_table_output_dir, mcl_inputs_tmp_dir):
    logger.info('Preparing jobs inputs for prepare_og_for_mcl...')

    job_index_to_ogs = {i: [] for i in range(max_parallel_jobs)}
    job_index_to_genes_count = {i: 0 for i in range(max_parallel_jobs)}

    ogs_genes_count = pd.read_csv(os.path.join(putative_orthologs_table_output_dir, 'num_of_genes_in_putative_sets.csv'))
    ogs_genes_count.sort_values(by='genes_count', ascending=False, inplace=True)
    for _, row in ogs_genes_count.iterrows():
        # Find the job with the smallest current genes count
        job_index_with_min_genes_count = min(job_index_to_genes_count, key=job_index_to_genes_count.get)
        # Assign the current OG to this job
        job_index_to_ogs[job_index_with_min_genes_count].append(row["OG_name"])
        # Update the job's genes count
        job_index_to_genes_count[job_index_with_min_genes_count] += row["genes_count"]

    job_inputs_dir = os.path.join(mcl_inputs_tmp_dir, 'jobs_inputs')
    os.makedirs(job_inputs_dir, exist_ok=True)
    for job_index, ogs in job_index_to_ogs.items():
        job_path = os.path.join(job_inputs_dir, f'{job_index}.txt')
        with open(job_path, 'w') as f:
            f.write('\n'.join(map(str, ogs)))

    logger.info('Done preparing jobs inputs for prepare_og_for_mcl.')


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    main(logger, MAX_PARALLEL_JOBS, SCRIPT_DIR, os.path.join(SCRIPT_DIR, 'mcl_inputs_tmp_dir'))
