import os
import math

from . import consts
from .pipeline_auxiliaries import wait_for_results, prepare_directories, fail, submit_mini_batch, submit_batch, execute
from .logic_auxiliaries import define_intervals, aggregate_mmseqs_scores
from .file_writer import write_to_file


def unify_clusters_mmseqs_hits(logger, times_logger, output_dir, tmp_dir, done_files_dir, error_file_path,
                               mmseqs_output_dir, run_optimized_mmseqs, queue_name, account_name):
    if run_optimized_mmseqs:
        rbhs_dir_name = '05bb_extract_rbh_hits'
        paraloogs_dir_name = '05bc_paralogs'
    else:
        rbhs_dir_name = '05ba_all_vs_all_analysis'
        paraloogs_dir_name = '05bc_mmseqs_paralogs'

    # 05_3. concat_clusters_info.py
    step_number = f'05_3'
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

    # 05_4. concat_clusters_info.py
    step_number = f'05_4'
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
                                       start_substep_number, account_name, queue_name):
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
                                 os.path.join(normalized_hits_output_dir, f"{strains_pair}.{step_name}"),
                                 scores_normalize_coefficients[strains_pair]
                                 ]
            all_cmds_params.append(single_cmd_params)

        for hits_file in os.listdir(paralogs_output_dir):
            if not hits_file.endswith('m8_filtered'):
                continue
            strains_pair = os.path.splitext(hits_file)[0]
            single_cmd_params = [os.path.join(paralogs_output_dir, hits_file),
                                 os.path.join(normalized_hits_output_dir, f"{strains_pair}.{step_name}"),
                                 scores_normalize_coefficients[strains_pair]
                                 ]
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

        number_of_hits_files = len(os.listdir(normalized_hits_output_dir))
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
