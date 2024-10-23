import os
import sys
import itertools
import json
from sys import argv
import argparse
import traceback

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, prepare_directories, submit_batch, \
    wait_for_results, submit_mini_batch, execute, get_job_times_logger
from auxiliaries.logic_auxiliaries import aggregate_mmseqs_scores, define_intervals
from auxiliaries import consts
from auxiliaries.file_writer import write_to_file


def full_orthogroups_infernece(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                               translated_orfs_dir, strains_names_path, queue_name,
                               account_name, identity_cutoff, coverage_cutoff, e_value_cutoff, num_of_times_script_called):
    with open(strains_names_path) as f:
        strains_names = f.read().rstrip().split('\n')
    n_jobs_per_step = consts.MAX_PARALLEL_JOBS // num_of_times_script_called

    # a.	mmseqs2_all_vs_all.py
    # Input: (1) 2 input paths for 2 (different) genome files (query and target), g1 and g2
    #        (2) an output file path (with a suffix as follows: i_vs_j.tsv. especially relevant for the wrapper).
    # Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
    # Precisely, for each gene x in g1, query x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
    # Can be parallelized on cluster
    step_number = f'{base_step_number}a'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_all_vs_all_analysis'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs2_all_vs_all.py')
    orthologs_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orthologs_scores_statistics_dir = os.path.join(orthologs_output_dir, 'scores_statistics')
    os.makedirs(orthologs_scores_statistics_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        orfs_files = [filename for filename in os.listdir(translated_orfs_dir) if filename.endswith('02d_translated_orfs')]
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
                                     error_file_path,
                                     f'--identity_cutoff {identity_cutoff / 100}',
                                     f'--coverage_cutoff {coverage_cutoff / 100}',
                                     f'--e_value_cutoff {e_value_cutoff}',
                                     ]

                all_cmds_params.append(single_cmd_params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                       num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                       job_name_suffix='rbh_analysis',
                                                       queue_name=queue_name,
                                                       account_name=account_name,
                                                       memory=consts.MMSEQS_REQUIRED_MEMORY_GB,
                                                       time_in_hours=72)

            wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                             num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # b. get_max_score_per_gene
    # Input: path to folder with all reciprocal hits files
    # Output: Dictionary with max score per gene
    # CANNOT be parallelized on cluster
    step_number = f'{base_step_number}b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    script_path = os.path.join(consts.SRC_DIR, 'steps/max_rbh_score.py')
    max_rbh_scores_step_name = f'{step_number}_max_rbh_scores'
    max_rbh_scores_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, max_rbh_scores_step_name)
    done_file_path = os.path.join(done_files_dir, f'{max_rbh_scores_step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Calculating max rbh scores per gene...')

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for strain_name in strains_names:
            single_cmd_params = [orthologs_output_dir, strain_name, max_rbh_scores_output_dir, max_rbh_scores_step_name,
                                 error_file_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='max_rbh_score',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, max_rbh_scores_step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # c. mmseqs2_paralogs.py
    # Input: max orthologs score for each gene
    # Output: paralogs in each genome
    step_number = f'{base_step_number}c'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mmseqs_paralogs'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs2_paralogs.py')
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
                                 error_file_path,
                                 f'--identity_cutoff {identity_cutoff / 100}',
                                 f'--coverage_cutoff {coverage_cutoff / 100}',
                                 f'--e_value_cutoff {e_value_cutoff}',
                                 ]

            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
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

    # d. normalize_scores
    step_number = f'{base_step_number}d'
    script_path = os.path.join(consts.SRC_DIR, 'steps/normalize_hits_scores.py')
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_normalize_scores'
    normalized_hits_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
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

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='hits_normalize',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # e. concatenate_all_hits
    # Input: path to folder with all hits files
    # Output: concatenated file of all hits files
    # CANNOT be parallelized on cluster
    step_number = f'{base_step_number}e'
    script_path = os.path.join(consts.SRC_DIR, 'steps/concatenate_hits.py')
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_concatenate_hits'
    concatenate_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    all_hits_file = os.path.join(concatenate_output_dir, 'concatenated_all_hits.txt')
    if not os.path.exists(done_file_path):
        logger.info('Concatenating hits...')

        number_of_hits_files = len(os.listdir(normalized_hits_output_dir))
        intervals = define_intervals(0, number_of_hits_files, n_jobs_per_step)

        concatenated_chunks_dir = os.path.join(concatenate_output_dir, 'temp_chunks')
        os.makedirs(concatenated_chunks_dir, exist_ok=True)
        all_cmds_params = [[normalized_hits_output_dir, start, end, concatenated_chunks_dir]
                           for (start, end) in intervals]

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='concatenate_hits',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        execute(logger, f'cat {concatenated_chunks_dir}/* >> {all_hits_file}', process_is_string=True)
        execute(logger, f'rm -rf {concatenated_chunks_dir}', process_is_string=True)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # f.	construct_putative_orthologs_table.py
    # Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
    # Output: updates the table with the info from the reciprocal hit file.
    # CANNOT be parallelized on cluster
    step_number = f'{base_step_number}f'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_putative_table'
    script_path = os.path.join(consts.SRC_DIR, 'steps/construct_putative_orthologs_table.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    putative_orthologs_table_path = os.path.join(pipeline_step_output_dir, 'putative_orthologs_table.txt')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing putative orthologs table...')
        job_name = os.path.split(script_path)[-1]
        params = [all_hits_file,
                  putative_orthologs_table_path]
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path, queue_name, account_name, job_name=job_name)
        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # g.   prepare_files_for_mcl.py
    # Input: (1) a path for a concatenated all reciprocal hits file (2) a path for a putative orthologs file (3) a path for an output folder
    # Output: an input file for MCL for each putative orthologs group
    # CANNOT be parallelized on cluster (if running on the concatenated file)
    step_number = f'{base_step_number}g'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mcl_input_files'
    script_path = os.path.join(consts.SRC_DIR, 'steps/prepare_files_for_mcl.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Preparing files for MCL...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        clusters_to_prepare_per_job = 100
        with open(os.path.join(os.path.split(putative_orthologs_table_path)[0], 'num_of_putative_sets.txt')) as f:
            num_of_putative_sets = int(f.read())
        if num_of_putative_sets == 0:
            error_msg = f'No putative ortholog groups were detected in your dataset. Please try to lower the ' \
                        f'similarity parameters (see Advanced Options in the submission page) and re-submit your job.'
            fail(logger, error_msg, error_file_path)
        for i in range(1, num_of_putative_sets + 1, clusters_to_prepare_per_job):
            first_mcl = str(i)
            last_mcl = str(min(i + clusters_to_prepare_per_job - 1,
                               num_of_putative_sets))  # 1-10, 11-20, etc... for clusters_to_prepare_per_job=10

            single_cmd_params = [all_hits_file,
                                 putative_orthologs_table_path,
                                 first_mcl,
                                 last_mcl,
                                 pipeline_step_output_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   # *times* the number of clusters_to_prepare_per_job above. 50 in total per batch!
                                                   job_name_suffix='mcl_preparation',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # h.	run_mcl.py
    # Input: (1) a path to an MCL input file (2) a path to MCL's output.
    # Output: MCL analysis.
    # Can be parallelized on cluster
    step_number = f'{base_step_number}h'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_mcl_analysis'
    script_path = os.path.join(consts.SRC_DIR, 'steps/run_mcl.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Executing MCL...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for putative_orthologs_group in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_group_prefix = os.path.splitext(putative_orthologs_group)[0]
            output_file_name = f'{putative_orthologs_group_prefix}.{step_name}'

            single_cmd_params = [f'"{os.path.join(previous_pipeline_step_output_dir, putative_orthologs_group)}"',
                                 f'"{os.path.join(pipeline_step_output_dir, output_file_name)}"']
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='mcl_execution',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # i.	verify_cluster.py
    # Input: (1) mcl analysis file (2) a path to which the file will be moved if relevant (3) optional: maximum number of clusters allowed [default=1]
    # Output: filter irrelevant clusters by moving the relevant to an output directory
    # Can be parallelized on cluster
    step_number = f'{base_step_number}i'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_verified_clusters'
    script_path = os.path.join(consts.SRC_DIR, 'steps/verify_cluster.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Verifying clusters...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for putative_orthologs_group in os.listdir(previous_pipeline_step_output_dir):
            single_cmd_params = [os.path.join(previous_pipeline_step_output_dir, putative_orthologs_group),
                                 pipeline_step_output_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='clusters_verification',
                                                   queue_name=queue_name,
                                                   account_name=account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    logger.info(f'A total of {len(os.listdir(previous_pipeline_step_output_dir))} clusters were analyzed. '
                f'{len(os.listdir(pipeline_step_output_dir))} clusters were produced.')

    # j.	construct_verified_orthologs_table.py
    # Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
    # Output: aggregates all the well-clustered OGs to the final table.
    step_number = f'{base_step_number}j'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_verified_table'
    script_path = os.path.join(consts.SRC_DIR, 'steps/construct_verified_orthologs_table.py')
    verified_orthologs_table_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    verified_orthologs_table_file_path = os.path.join(verified_orthologs_table_dir_path, 'verified_orthologs_table.csv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing verified orthologs table...')
        params = [putative_orthologs_table_path,
                  previous_pipeline_step_output_dir,
                  verified_orthologs_table_file_path]

        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path, queue_name, account_name,
                          job_name='verified_ortholog_groups')
        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path,
                         error_message='No ortholog groups were detected in your dataset. Please try to lower '
                                       'the similarity parameters (see Advanced Options in the submission page) '
                                       'and re-submit your job.')

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('step_number', help='')
    parser.add_argument('output_dir', help='')
    parser.add_argument('tmp_dir', help='')
    parser.add_argument('done_files_dir', help='')
    parser.add_argument('translated_orfs_dir', help='')
    parser.add_argument('strains_names_path', help='')
    parser.add_argument('queue_name', help='')
    parser.add_argument('account_name', help='')
    parser.add_argument('identity_cutoff', help='', type=float)
    parser.add_argument('coverage_cutoff', help='', type=float)
    parser.add_argument('e_value_cutoff', help='', type=float)
    parser.add_argument('num_of_times_script_called', help='', type=int)
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir)
    times_logger = get_job_times_logger(args.logs_dir)

    logger.info(script_run_message)
    try:
        full_orthogroups_infernece(logger, times_logger, args.step_number, args.error_file_path, args.output_dir,
                                   args.tmp_dir, args.done_files_dir, args.translated_orfs_dir, args.strains_names_path,
                                   args.queue_name, args.account_name,
                                   args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff, args.num_of_times_script_called)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
