import os
import math
from collections import defaultdict
import dask.dataframe as dd
import itertools

from . import consts
from .pipeline_auxiliaries import wait_for_results, prepare_directories, fail, submit_mini_batch, submit_batch, execute
from .logic_auxiliaries import aggregate_mmseqs_scores, define_intervals, add_score_column_to_mmseqs_output, get_directory_size_in_gb
from .file_writer import write_to_file


def run_mmseqs_and_extract_hits(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                               translated_orfs_dir, all_proteins_path, strains_names_path, queue_name,
                               account_name, identity_cutoff, coverage_cutoff, e_value_cutoff, max_parallel_jobs,
                               run_optimized_mmseqs, use_parquet, verbose):
    with open(strains_names_path) as f:
        strains_names = f.read().rstrip().split('\n')

    if run_optimized_mmseqs:
        orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir =\
            run_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir,
                           done_files_dir, all_proteins_path, strains_names, queue_name, account_name,
                           identity_cutoff, coverage_cutoff,  e_value_cutoff, max_parallel_jobs, use_parquet, verbose)
    else:
        orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir =\
            run_non_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir,
                               done_files_dir, translated_orfs_dir, strains_names, queue_name, account_name,
                               identity_cutoff, coverage_cutoff, e_value_cutoff, max_parallel_jobs, use_parquet)

    return orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir


def unify_clusters_mmseqs_hits(logger, times_logger, output_dir, tmp_dir, done_files_dir, error_file_path,
                               mmseqs_output_dir, run_optimized_mmseqs, queue_name, account_name, base_step_number,
                               start_substep_number):
    if run_optimized_mmseqs:
        rbhs_dir_name = '05_1_2_extract_rbh_hits'
        paraloogs_dir_name = '05_1_3_paralogs'
    else:
        rbhs_dir_name = '05_1_1_all_vs_all_analysis'
        paraloogs_dir_name = '05_1_3_mmseqs_paralogs'

    #  concat_clusters_info.py
    step_number = f'{base_step_number}_{start_substep_number}'
    script_path = os.path.join(consts.SRC_DIR, 'steps/concat_clusters_info.py')
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_concat_clusters_rbhs'
    orthologs_output_dir, orthologs_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orthologs_scores_statistics_dir = os.path.join(orthologs_output_dir, 'scores_statistics')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        if len(os.listdir(mmseqs_output_dir)) == 1:
            cmd = f'cp -r {os.path.join(mmseqs_output_dir, "0", rbhs_dir_name)} {orthologs_output_dir}'
            execute(logger, cmd, process_is_string=True)
        else:
            os.makedirs(orthologs_scores_statistics_dir, exist_ok=True)
            all_cmds_params = []
            hits_files_processed = set()
            for cluster_dir in os.listdir(mmseqs_output_dir):
                for hits_file in os.listdir(os.path.join(mmseqs_output_dir, cluster_dir, rbhs_dir_name)):
                    if not hits_file.endswith('m8') or hits_file in hits_files_processed:
                        continue
                    single_cmd_params = [mmseqs_output_dir, rbhs_dir_name, hits_file, orthologs_output_dir]
                    all_cmds_params.append(single_cmd_params)
                    hits_files_processed.add(hits_file)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, orthologs_tmp_dir, error_file_path,
                                                       num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                       job_name_suffix='concat_clusters_rbhs',
                                                       queue_name=queue_name,
                                                       account_name=account_name)

            wait_for_results(logger, times_logger, step_name, orthologs_tmp_dir, num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # concat_clusters_info.py
    step_number = f'{base_step_number}_{start_substep_number + 1}'
    script_path = os.path.join(consts.SRC_DIR, 'steps/concat_clusters_info.py')
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_concat_clusters_paralogs'
    paralogs_output_dir, paralogs_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    paralogs_scores_statistics_dir = os.path.join(paralogs_output_dir, 'scores_statistics')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        if len(os.listdir(mmseqs_output_dir)) == 1:
            cmd = f'cp -r {os.path.join(mmseqs_output_dir, "0", paraloogs_dir_name)} {paralogs_output_dir}'
            execute(logger, cmd, process_is_string=True)
        else:
            os.makedirs(paralogs_scores_statistics_dir, exist_ok=True)
            all_cmds_params = []
            hits_files_processed = set()
            for cluster_dir in os.listdir(mmseqs_output_dir):
                for hits_file in os.listdir(os.path.join(mmseqs_output_dir, cluster_dir, paraloogs_dir_name)):
                    if not hits_file.endswith('m8_filtered') or hits_file in hits_files_processed:
                        continue
                    single_cmd_params = [mmseqs_output_dir, paraloogs_dir_name, hits_file, paralogs_output_dir]
                    all_cmds_params.append(single_cmd_params)
                    hits_files_processed.add(hits_file)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, paralogs_tmp_dir, error_file_path,
                                                       num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                       job_name_suffix='concat_clusters_rbhs',
                                                       queue_name=queue_name,
                                                       account_name=account_name)

            wait_for_results(logger, times_logger, step_name, paralogs_tmp_dir, num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthologs_output_dir, orthologs_scores_statistics_dir, paralogs_output_dir, paralogs_scores_statistics_dir


def cluster_mmseqs_hits_to_orthogroups(logger, times_logger, error_file_path, output_dir, tmp_dir, done_files_dir,
                                       orthologs_output_dir, orthologs_scores_statistics_dir, paralogs_output_dir,
                                       paralogs_scores_statistics_dir, max_parallel_jobs, base_step_number,
                                       start_substep_number, account_name, queue_name, use_parquet):
    # normalize_scores
    step_number = f'{base_step_number}_{start_substep_number}'
    script_path = os.path.join(consts.SRC_DIR, 'steps/normalize_hits_scores.py')
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_normalize_scores'
    normalized_hits_output_dir, normalized_hits_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        scores_statistics_file = os.path.join(normalized_hits_output_dir, 'scores_stats.json')
        scores_normalize_coefficients = aggregate_mmseqs_scores(orthologs_scores_statistics_dir, paralogs_scores_statistics_dir, scores_statistics_file)

        all_cmds_params = []
        for hits_file in os.listdir(orthologs_output_dir):
            if not hits_file.endswith('m8'):
                continue
            strains_pair = os.path.splitext(hits_file)[0]
            single_cmd_params = [os.path.join(orthologs_output_dir, hits_file),
                                 os.path.join(normalized_hits_output_dir, f"{strains_pair}.m8"),
                                 scores_normalize_coefficients[strains_pair]
                                 ]
            if use_parquet:
                single_cmd_params.append('--use_parquet')
            all_cmds_params.append(single_cmd_params)

        for hits_file in os.listdir(paralogs_output_dir):
            if not hits_file.endswith('m8_filtered'):
                continue
            strains_pair = os.path.splitext(hits_file)[0]
            single_cmd_params = [os.path.join(paralogs_output_dir, hits_file),
                                 os.path.join(normalized_hits_output_dir, f"{strains_pair}.m8"),
                                 scores_normalize_coefficients[strains_pair]
                                 ]
            if use_parquet:
                single_cmd_params.append('--use_parquet')
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, normalized_hits_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // max_parallel_jobs),
                                                   job_name_suffix='hits_normalize',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, normalized_hits_tmp_dir, num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # concatenate_all_hits
    # Input: path to folder with all hits files
    # Output: concatenated file of all hits files
    # CANNOT be parallelized on cluster
    step_number = f'{base_step_number}_{start_substep_number + 1}'
    script_path = os.path.join(consts.SRC_DIR, 'steps/concatenate_hits.py')
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_concatenate_hits'
    concatenate_output_dir, concatenate_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    all_hits_file = os.path.join(concatenate_output_dir, 'concatenated_all_hits.txt')
    if not os.path.exists(done_file_path):
        logger.info('Concatenating hits...')

        number_of_hits_files = sum(1 for file in os.listdir(normalized_hits_output_dir) if file.endswith('.m8'))
        intervals = define_intervals(0, number_of_hits_files, max_parallel_jobs)

        concatenated_chunks_dir = os.path.join(concatenate_output_dir, 'temp_chunks')
        os.makedirs(concatenated_chunks_dir, exist_ok=True)
        all_cmds_params = [[normalized_hits_output_dir, start, end, concatenated_chunks_dir]
                           for (start, end) in intervals]

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, concatenate_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // max_parallel_jobs),
                                                   job_name_suffix='concatenate_hits',
                                                    queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, concatenate_tmp_dir,
                         num_of_batches, error_file_path)

        execute(logger, f'cat {concatenated_chunks_dir}/* >> {all_hits_file}', process_is_string=True)
        execute(logger, f'rm -rf {concatenated_chunks_dir}', process_is_string=True)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # construct_putative_orthologs_table.py
    # Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
    # Output: updates the table with the info from the reciprocal hit file.
    # CANNOT be parallelized on cluster
    step_number = f'{base_step_number}_{start_substep_number + 2}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_putative_table'
    script_path = os.path.join(consts.SRC_DIR, 'steps/construct_putative_orthologs_table.py')
    putative_orthologs_table_output_dir, putative_orthologs_table_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    putative_orthologs_table_path = os.path.join(putative_orthologs_table_output_dir, 'putative_orthologs_table.txt')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing putative orthologs table...')
        job_name = os.path.split(script_path)[-1]
        params = [all_hits_file,
                  putative_orthologs_table_path]
        submit_mini_batch(logger, script_path, [params], putative_orthologs_table_tmp_dir, error_file_path,
                          queue_name, account_name, job_name=job_name)
        wait_for_results(logger, times_logger, step_name, putative_orthologs_table_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # prepare_files_for_mcl.py
    # Input: (1) a path for a concatenated all reciprocal hits file (2) a path for a putative orthologs file (3) a path for an output folder
    # Output: an input file for MCL for each putative orthologs group
    # CANNOT be parallelized on cluster (if running on the concatenated file)
    step_number = f'{base_step_number}_{start_substep_number + 3}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mcl_input_files'
    script_path = os.path.join(consts.SRC_DIR, 'steps/prepare_files_for_mcl.py')
    mcl_inputs_dir, mcl_inputs_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Preparing files for MCL...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs

        with open(os.path.join(os.path.split(putative_orthologs_table_path)[0], 'num_of_putative_sets.txt')) as f:
            num_of_putative_sets = int(f.read())
        if num_of_putative_sets == 0:
            error_msg = f'No putative ortholog groups were detected in your dataset. Please try to lower the ' \
                        f'similarity parameters (see Advanced Options in the submission page) and re-submit your job.'
            fail(logger, error_msg, error_file_path)
        clusters_to_prepare_per_job = math.ceil(num_of_putative_sets / max_parallel_jobs)
        for i in range(1, num_of_putative_sets + 1, clusters_to_prepare_per_job):
            first_mcl = str(i)
            last_mcl = str(min(i + clusters_to_prepare_per_job - 1,
                               num_of_putative_sets))  # 1-10, 11-20, etc... for clusters_to_prepare_per_job=10

            single_cmd_params = [all_hits_file,
                                 putative_orthologs_table_path,
                                 first_mcl,
                                 last_mcl,
                                 mcl_inputs_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, mcl_inputs_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // max_parallel_jobs),
                                                   # *times* the number of clusters_to_prepare_per_job above. 50 in total per batch!
                                                   job_name_suffix='mcl_preparation',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, mcl_inputs_tmp_dir, num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # run_mcl.py
    # Input: (1) a path to an MCL input file (2) a path to MCL's output.
    # Output: MCL analysis.
    # Can be parallelized on cluster
    step_number = f'{base_step_number}_{start_substep_number + 4}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mcl_analysis'
    script_path = os.path.join(consts.SRC_DIR, 'steps/run_mcl.py')
    mcl_outputs_dir, mcl_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Executing MCL...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for putative_orthologs_group in os.listdir(mcl_inputs_dir):
            putative_orthologs_group_prefix = os.path.splitext(putative_orthologs_group)[0]
            output_file_name = f'{putative_orthologs_group_prefix}.{step_name}'

            single_cmd_params = [f'"{os.path.join(mcl_inputs_dir, putative_orthologs_group)}"',
                                 f'"{os.path.join(mcl_outputs_dir, output_file_name)}"']
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, mcl_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // max_parallel_jobs),
                                                   job_name_suffix='mcl_execution',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, mcl_tmp_dir,
                         num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # verify_cluster.py
    # Input: (1) mcl analysis file (2) a path to which the file will be moved if relevant (3) optional: maximum number of clusters allowed [default=1]
    # Output: filter irrelevant clusters by moving the relevant to an output directory
    # Can be parallelized on cluster
    step_number = f'{base_step_number}_{start_substep_number + 5}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_verified_clusters'
    script_path = os.path.join(consts.SRC_DIR, 'steps/verify_cluster.py')
    verified_clusters_output_dir, verified_clusters_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Verifying clusters...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for putative_orthologs_group in os.listdir(mcl_outputs_dir):
            single_cmd_params = [os.path.join(mcl_outputs_dir, putative_orthologs_group),
                                 verified_clusters_output_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, verified_clusters_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // max_parallel_jobs),
                                                   job_name_suffix='clusters_verification',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, verified_clusters_tmp_dir, num_of_batches, error_file_path)

        logger.info(f'A total of {mcl_inputs_dir} clusters were analyzed. '
                    f'{len(os.listdir(verified_clusters_output_dir))} clusters were produced.')

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # construct_verified_orthologs_table.py
    # Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
    # Output: aggregates all the well-clustered OGs to the final table.
    step_number = f'{base_step_number}_{start_substep_number + 6}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_verified_table'
    script_path = os.path.join(consts.SRC_DIR, 'steps/construct_verified_orthologs_table.py')
    verified_orthologs_table_dir_path, verified_orthologs_table_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orthogroups_file_path = os.path.join(verified_orthologs_table_dir_path, 'orthogroups.csv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing verified orthologs table...')
        params = [putative_orthologs_table_path,
                  verified_clusters_output_dir,
                  orthogroups_file_path]

        submit_mini_batch(logger, script_path, [params], verified_orthologs_table_tmp_dir, error_file_path,
                          queue_name, account_name, job_name='verified_ortholog_groups')
        wait_for_results(logger, times_logger, step_name, verified_orthologs_table_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path,
                         error_message='No ortholog groups were detected in your dataset. Please try to lower '
                                       'the similarity parameters (see Advanced Options in the submission page) '
                                       'and re-submit your job.')

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthogroups_file_path


def run_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                       all_proteins_path, strains_names, queue_name, account_name, identity_cutoff, coverage_cutoff,
                       e_value_cutoff, n_jobs_per_step, use_parquet, verbose):
    # 1.	mmseqs2_all_vs_all.py
    step_number = f'{base_step_number}_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_all_vs_all_analysis'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs_v2/mmseqs2_all_vs_all.py')
    all_vs_all_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    m8_output_path = os.path.join(all_vs_all_output_dir, 'all_vs_all_reduced_columns.csv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        m8_raw_output_path = os.path.join(all_vs_all_output_dir, 'all_vs_all_raw.m8')
        params = [all_proteins_path,
                  all_vs_all_output_dir,
                  m8_raw_output_path,
                  f'--identity_cutoff {identity_cutoff / 100}',
                  f'--coverage_cutoff {coverage_cutoff / 100}',
                  f'--e_value_cutoff {e_value_cutoff}',
                  f'--number_of_genomes {len(strains_names)}',
                  f'--cpus {consts.MMSEQS_NUM_OF_CORES}']

        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path, queue_name,
                          account_name, job_name='mmseqs', num_of_cpus=consts.MMSEQS_NUM_OF_CORES,
                          memory=consts.MMSEQS_BIG_DATASET_REQUIRED_MEMORY_GB, time_in_hours=72)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, 1, error_file_path)

        logger.info(f"Starting to read {m8_raw_output_path} and create a reduced version of it...")
        m8_df = dd.read_csv(m8_raw_output_path, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER, dtype=consts.MMSEQS_OUTPUT_COLUMNS_TYPES)
        m8_df = m8_df[m8_df['query'] != m8_df['target']]
        add_score_column_to_mmseqs_output(m8_df)
        m8_df['query_genome'] = m8_df['query'].str.split(':').str[0]
        m8_df['target_genome'] = m8_df['target'].str.split(':').str[0]
        m8_df = m8_df[['query', 'query_genome', 'target', 'target_genome', 'score']]
        m8_df.to_parquet(m8_output_path)  # Here I always use parquet since it's a huge file
        logger.info(f"{m8_output_path} was created successfully.")

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    m8_output_size_in_gb = get_directory_size_in_gb(m8_output_path)
    m8_output_parsing_memory = str(max(int(consts.DEFAULT_MEMORY_PER_JOB_GB), math.ceil(m8_output_size_in_gb * 10)))

    # 2.	extract_rbh_hits.py
    step_number = f'{base_step_number}_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_extract_rbh_hits'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs_v2/extract_rbh_hits.py')
    orthologs_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orthologs_scores_statistics_dir = os.path.join(orthologs_output_dir, 'scores_statistics')
    os.makedirs(orthologs_scores_statistics_dir, exist_ok=True)
    max_rbh_scores_parts_output_dir = os.path.join(orthologs_output_dir, 'max_rbh_scores_parts')
    os.makedirs(max_rbh_scores_parts_output_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        if len(strains_names) >= 2:
            rbh_inputs_dir = os.path.join(pipeline_step_tmp_dir, 'rbh_inputs_dir')
            os.makedirs(rbh_inputs_dir, exist_ok=True)
            job_index_to_pairs = defaultdict(list)
            for i, (genome1, genome2) in enumerate(itertools.combinations(strains_names, 2)):
                job_index = i % n_jobs_per_step
                job_index_to_pairs[job_index].append((genome1, genome2))

            for job_index, pairs in job_index_to_pairs.items():
                job_input_path = os.path.join(rbh_inputs_dir, f'job_{job_index}_pairs.txt')
                pairs_text = '\n'.join([f'{genome1}\t{genome2}' for genome1, genome2 in pairs])
                with open(job_input_path, 'w') as f:
                    f.write(pairs_text)

            all_cmds_params = []
            for rbh_input_file in os.listdir(rbh_inputs_dir):
                rbh_input_path = os.path.join(rbh_inputs_dir, rbh_input_file)
                params = [m8_output_path, rbh_input_path, orthologs_output_dir,
                          orthologs_scores_statistics_dir, max_rbh_scores_parts_output_dir]
                if use_parquet:
                    params.append('--use_parquet')
                if verbose:
                    params.append('--verbose')
                all_cmds_params.append(params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                       error_file_path,
                                                       num_of_cmds_per_job=1, # = max(1, len(all_cmds_params) // n_jobs_per_step),
                                                       job_name_suffix='rbh_analysis',
                                                       queue_name=queue_name,
                                                       account_name=account_name,
                                                       memory=m8_output_parsing_memory)

            wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                             num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 3. extract_paralogs.py
    step_number = f'{base_step_number}_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_paralogs'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs_v2/extract_paralogs.py')
    paralogs_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    paralogs_scores_statistics_dir = os.path.join(paralogs_output_dir, 'scores_statistics')
    os.makedirs(paralogs_scores_statistics_dir, exist_ok=True)
    max_rbh_scores_unified_dir = os.path.join(paralogs_output_dir, 'max_rbh_scores_unified')
    os.makedirs(max_rbh_scores_unified_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Searching for paralogs in each genome')

        paralogs_inputs_dir = os.path.join(pipeline_step_tmp_dir, 'paralogs_inputs_dir')
        os.makedirs(paralogs_inputs_dir, exist_ok=True)
        job_index_to_genome = defaultdict(list)
        for i, genome in enumerate(strains_names):
            job_index = i % n_jobs_per_step
            job_index_to_genome[job_index].append(genome)

        for job_index, genomes in job_index_to_genome.items():
            job_input_path = os.path.join(paralogs_inputs_dir, f'job_{job_index}_genomes.txt')
            genomes_text = '\n'.join(genomes)
            with open(job_input_path, 'w') as f:
                f.write(genomes_text)

        all_cmds_params = []
        for paralogs_input_file in os.listdir(paralogs_inputs_dir):
            paralogs_input_path = os.path.join(paralogs_inputs_dir, paralogs_input_file)
            single_cmd_params = [m8_output_path, paralogs_input_path, max_rbh_scores_parts_output_dir,
                                 paralogs_output_dir, max_rbh_scores_unified_dir, paralogs_scores_statistics_dir]
            if use_parquet:
                single_cmd_params.append('--use_parquet')
            if verbose:
                single_cmd_params.append('--verbose')
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   error_file_path,
                                                   num_of_cmds_per_job=1, # = max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='paralogs_analysis',
                                                   queue_name=queue_name,
                                                   account_name=account_name,
                                                   memory=m8_output_parsing_memory)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir


def run_non_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                           translated_orfs_dir, strains_names, queue_name, account_name, identity_cutoff,
                           coverage_cutoff, e_value_cutoff, n_jobs_per_step, use_parquet):
    # 1.	mmseqs2_all_vs_all.py
    # Input: (1) 2 input paths for 2 (different) genome files (query and target), g1 and g2
    #        (2) an output file path (with a suffix as follows: i_vs_j.tsv. especially relevant for the wrapper).
    # Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
    # Precisely, for each gene x in g1, query x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
    # Can be parallelized on cluster
    step_number = f'{base_step_number}_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_all_vs_all_analysis'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs_v1/mmseqs2_all_vs_all.py')
    orthologs_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orthologs_scores_statistics_dir = os.path.join(orthologs_output_dir, 'scores_statistics')
    os.makedirs(orthologs_scores_statistics_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        orfs_files = [filename for filename in os.listdir(translated_orfs_dir) if filename.endswith('.faa')]
        if len(orfs_files) >= 2:
            logger.info(f'Querying all VS all (using mmseqs2)...')
            all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
            for fasta1, fasta2 in itertools.combinations(orfs_files, 2):
                strain1_name = os.path.splitext(fasta1)[0]
                strain2_name = os.path.splitext(fasta2)[0]

                output_file_name = f'{strain1_name}_vs_{strain2_name}.m8'
                output_file_path = os.path.join(orthologs_output_dir, output_file_name)

                single_cmd_params = [os.path.join(translated_orfs_dir, fasta1),
                                     os.path.join(translated_orfs_dir, fasta2),
                                     output_file_path,
                                     orthologs_scores_statistics_dir,
                                     f'--identity_cutoff {identity_cutoff / 100}',
                                     f'--coverage_cutoff {coverage_cutoff / 100}',
                                     f'--e_value_cutoff {e_value_cutoff}',
                                     ]
                if use_parquet:
                    single_cmd_params.append('--use_parquet')
                all_cmds_params.append(single_cmd_params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                       error_file_path,
                                                       num_of_cmds_per_job=max(1,
                                                                               len(all_cmds_params) // n_jobs_per_step),
                                                       job_name_suffix='rbh_analysis',
                                                       queue_name=queue_name,
                                                       account_name=account_name,
                                                       time_in_hours=72)

            wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                             num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 2. get_max_score_per_gene
    # Input: path to folder with all reciprocal hits files
    # Output: Dictionary with max score per gene
    # CANNOT be parallelized on cluster
    step_number = f'{base_step_number}_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs_v1/max_rbh_score.py')
    max_rbh_scores_step_name = f'{step_number}_max_rbh_scores'
    max_rbh_scores_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir,
                                                                           max_rbh_scores_step_name)
    done_file_path = os.path.join(done_files_dir, f'{max_rbh_scores_step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Calculating max rbh scores per gene...')

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for strain_name in strains_names:
            single_cmd_params = [orthologs_output_dir, strain_name, max_rbh_scores_output_dir, max_rbh_scores_step_name]
            if use_parquet:
                single_cmd_params.append('--use_parquet')
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='max_rbh_score',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, max_rbh_scores_step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 3. mmseqs2_paralogs.py
    # Input: max orthologs score for each gene
    # Output: paralogs in each genome
    step_number = f'{base_step_number}_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mmseqs_paralogs'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs_v1/mmseqs2_paralogs.py')
    paralogs_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    paralogs_scores_statistics_dir = os.path.join(paralogs_output_dir, 'scores_statistics')
    os.makedirs(paralogs_scores_statistics_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Searching for paralogs in each genome')
        all_cmds_params = []
        for fasta_file in os.listdir(translated_orfs_dir):
            strain_name = os.path.splitext(fasta_file)[0]

            logger.info(f'Querying {fasta_file} for paralogs')
            output_file_name = f'{strain_name}_vs_{strain_name}.m8'
            output_file_path = os.path.join(paralogs_output_dir, output_file_name)
            max_scores_file_path = os.path.join(max_rbh_scores_output_dir, f"{strain_name}.{max_rbh_scores_step_name}")

            single_cmd_params = [os.path.join(translated_orfs_dir, fasta_file),
                                 output_file_path,
                                 max_scores_file_path,
                                 paralogs_scores_statistics_dir,
                                 f'--identity_cutoff {identity_cutoff / 100}',
                                 f'--coverage_cutoff {coverage_cutoff / 100}',
                                 f'--e_value_cutoff {e_value_cutoff}',
                                 ]
            if use_parquet:
                single_cmd_params.append('--use_parquet')
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='paralogs_analysis',
                                                   queue_name=queue_name,
                                                   account_name=account_name,
                                                   memory=consts.MMSEQS_REQUIRED_MEMORY_GB)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir
