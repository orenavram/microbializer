import os
import math
from collections import defaultdict
import dask.dataframe as dd
import pandas as pd
import itertools
import shutil
import time
from datetime import timedelta

from . import consts
from .pipeline_auxiliaries import wait_for_results, prepare_directories, submit_mini_batch, submit_batch, write_done_file
from .logic_auxiliaries import aggregate_mmseqs_scores, add_score_column_to_mmseqs_output, \
    get_directory_size_in_gb, combine_orphan_genes_stats, split_ogs_to_jobs_inputs_files_by_og_sizes


def infer_orthogroups(logger, times_logger, config, infer_orthogroups_config, base_step_number):

    orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir = \
        run_mmseqs_and_extract_hits(logger, times_logger, config, infer_orthogroups_config, base_step_number)

    final_orthogroups_dir_path, orphan_genes_dir, final_substep_number = \
        cluster_mmseqs_hits_to_orthogroups(logger, times_logger, config, infer_orthogroups_config, orthologs_output_dir,
                                           orthologs_scores_statistics_dir, paralogs_output_dir,
                                           paralogs_scores_statistics_dir, base_step_number, 4)

    return final_orthogroups_dir_path, orphan_genes_dir, final_substep_number


def run_mmseqs_and_extract_hits(logger, times_logger, config, infer_orthogroups_config, base_step_number):
    if config.run_optimized_mmseqs:
        orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir =\
            run_unified_mmseqs(logger, times_logger, config, infer_orthogroups_config, base_step_number)
    else:
        orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir = \
            run_non_unified_mmseqs_with_dbs(logger, times_logger, config, infer_orthogroups_config, base_step_number)

    return orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir


def cluster_mmseqs_hits_to_orthogroups(logger, times_logger, config, infer_orthogroups_config, orthologs_output_dir,
                                       orthologs_scores_statistics_dir, paralogs_output_dir,
                                       paralogs_scores_statistics_dir, base_step_number, start_substep_number):
    # normalize_scores
    step_number = f'{base_step_number}_{start_substep_number}'
    script_path = consts.SRC_DIR / 'steps' / 'normalize_hits_scores.py'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_normalize_scores'
    normalized_hits_output_dir, normalized_hits_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Start normalizing hits scores...')
        start_time = time.time()

        scores_statistics_file = normalized_hits_output_dir / 'scores_stats.json'
        scores_normalize_coefficients = aggregate_mmseqs_scores(
            orthologs_scores_statistics_dir, paralogs_scores_statistics_dir, scores_statistics_file,
            infer_orthogroups_config.skip_paralogs)

        job_index_to_hits_files_and_coefficients = defaultdict(list)

        for i, hits_file in enumerate(orthologs_output_dir.glob('*.m8')):
            job_index = i % infer_orthogroups_config.max_parallel_jobs
            strains_pair = hits_file.stem
            job_index_to_hits_files_and_coefficients[job_index].append((str(hits_file),
                                                                        scores_normalize_coefficients[strains_pair]))

        if not infer_orthogroups_config.skip_paralogs:
            for i, hits_file in enumerate(paralogs_output_dir.glob('*.m8_filtered')):
                job_index = i % infer_orthogroups_config.max_parallel_jobs
                strains_pair = hits_file.stem
                job_index_to_hits_files_and_coefficients[job_index].append((str(hits_file),
                                                                            scores_normalize_coefficients[strains_pair]))

        all_cmds_params = []
        jobs_inputs_dir = normalized_hits_tmp_dir / 'jobs_inputs'
        os.makedirs(jobs_inputs_dir, exist_ok=True)
        for job_index, job_input_info in job_index_to_hits_files_and_coefficients.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                for hits_file, scores_normalize_coefficient in job_input_info:
                    f.write(f'{hits_file}\t{scores_normalize_coefficient}\n')

            single_cmd_params = [job_input_path, normalized_hits_output_dir]
            if config.use_parquet:
                single_cmd_params.append('--use_parquet')
            all_cmds_params.append(single_cmd_params)

        step_pre_processing_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} pre-processing took {step_pre_processing_time}.')

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, normalized_hits_tmp_dir,
                                      'hits_normalize')

        wait_for_results(logger, times_logger, step_name, normalized_hits_tmp_dir, num_of_batches, config.error_file_path)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


    # construct_putative_orthologs_table.py
    # Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
    # Output: updates the table with the info from the reciprocal hit file.
    # CANNOT be parallelized on cluster
    step_number = f'{base_step_number}_{start_substep_number + 1}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_putative_table'
    script_path = consts.SRC_DIR / 'steps' / 'construct_putative_orthologs_table.py'
    putative_orthologs_table_output_dir, putative_orthologs_table_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    putative_orthologs_table_path = putative_orthologs_table_output_dir / 'putative_orthologs_table.csv'
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Constructing putative orthologs table...')
        params = [normalized_hits_output_dir,
                  putative_orthologs_table_path]
        submit_mini_batch(logger, config, script_path, [params], putative_orthologs_table_tmp_dir,
                          os.path.split(script_path)[-1])
        wait_for_results(logger, times_logger, step_name, putative_orthologs_table_tmp_dir, 1,
                         config.error_file_path)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    putative_orthologs_table_df = pd.read_csv(putative_orthologs_table_path)

    # prepare_og_for_mcl.py
    step_number = f'{base_step_number}_{start_substep_number + 2}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mcl_input_files'
    script_path = consts.SRC_DIR / 'steps' / 'prepare_og_for_mcl.py'
    mcl_inputs_dir, mcl_inputs_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Preparing jobs inputs for prepare_og_for_mcl...')
        start_time = time.time()

        job_paths = split_ogs_to_jobs_inputs_files_by_og_sizes(putative_orthologs_table_df, mcl_inputs_tmp_dir,
                                                               infer_orthogroups_config.max_parallel_jobs)
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_path in job_paths:
            single_cmd_params = [normalized_hits_output_dir, putative_orthologs_table_path, job_path, mcl_inputs_dir]
            all_cmds_params.append(single_cmd_params)

        logger.info('Done preparing jobs inputs for prepare_og_for_mcl.')
        step_pre_processing_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} pre-processing time took {step_pre_processing_time}.')

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, mcl_inputs_tmp_dir,
                                      'mcl_preparation')

        wait_for_results(logger, times_logger, step_name, mcl_inputs_tmp_dir, num_of_batches, config.error_file_path)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


    # run_mcl.py
    # Input: (1) a path to an MCL input file (2) a path to MCL's output.
    # Output: MCL analysis.
    # Can be parallelized on cluster
    step_number = f'{base_step_number}_{start_substep_number + 3}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mcl_analysis'
    script_path = consts.SRC_DIR / 'steps' / 'run_mcl.py'
    mcl_step_outputs_dir, mcl_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    mcl_outputs_dir = mcl_step_outputs_dir / 'mcl_outputs'
    verified_clusters_output_dir = mcl_step_outputs_dir / 'verified_clusters'
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Executing MCL...')
        start_time = time.time()

        os.makedirs(mcl_outputs_dir, exist_ok=True)
        os.makedirs(verified_clusters_output_dir, exist_ok=True)

        job_paths = split_ogs_to_jobs_inputs_files_by_og_sizes(putative_orthologs_table_df, mcl_tmp_dir,
                                                               infer_orthogroups_config.max_parallel_jobs)
        all_cmds_params = []
        for job_path in job_paths:
            single_cmd_params = [mcl_inputs_dir, job_path, mcl_outputs_dir, verified_clusters_output_dir]
            all_cmds_params.append(single_cmd_params)

        logger.info('Done preparing jobs inputs for mcl_analysis.')
        step_pre_processing_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} pre-processing time took {step_pre_processing_time}.')

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, mcl_tmp_dir,
                                      'mcl_execution')

        wait_for_results(logger, times_logger, step_name, mcl_tmp_dir,
                         num_of_batches, config.error_file_path)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


    # construct_verified_orthologs_table.py
    # Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
    # Output: aggregates all the well-clustered OGs to the final table.
    step_number = f'{base_step_number}_{start_substep_number + 4}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_verified_table'
    script_path = consts.SRC_DIR / 'steps' 'construct_verified_orthologs_table.py'
    verified_orthologs_table_dir_path, verified_orthologs_table_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    orthogroups_file_path = verified_orthologs_table_dir_path / 'orthogroups.csv'
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not os.path.exists(done_file_path):
        logger.info('Constructing verified orthologs table...')
        params = [putative_orthologs_table_path,
                  verified_clusters_output_dir,
                  orthogroups_file_path]

        submit_mini_batch(logger, config, script_path, [params], verified_orthologs_table_tmp_dir,
                          'verified_ortholog_groups')
        wait_for_results(logger, times_logger, step_name, verified_orthologs_table_tmp_dir,
                         num_of_expected_results=1, error_file_path=config.error_file_path,
                         error_message='No ortholog groups were detected in your dataset. Please try to lower '
                                       'the similarity parameters (see Advanced Options in the submission page) '
                                       'and re-submit your job.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


    # extract_orphan_genes.py
    step_number = f'{base_step_number}_{start_substep_number + 5}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orphan_genes'
    script_path = consts.SRC_DIR / 'steps' / 'extract_orphan_genes.py'
    orphan_genes_dir, orphans_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    orphan_genes_internal_dir = orphan_genes_dir / 'orphans_lists_per_genome'
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Extracting orphan genes...')
        os.makedirs(orphan_genes_internal_dir, exist_ok=True)

        job_index_to_fasta_files = defaultdict(list)

        for i, fasta_file in enumerate(infer_orthogroups_config.translated_orfs_dir.iterdir()):
            job_index = i % infer_orthogroups_config.max_parallel_jobs
            job_index_to_fasta_files[job_index].append(str(fasta_file))

        jobs_inputs_dir = orphans_tmp_dir / 'job_inputs'
        os.makedirs(jobs_inputs_dir, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_fasta_files in job_index_to_fasta_files.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_fasta_files))

            single_cmd_params = [job_input_path, orthogroups_file_path, orphan_genes_internal_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, orphans_tmp_dir,
                                      'extract_orphans')

        wait_for_results(logger, times_logger, step_name, orphans_tmp_dir,
                         num_of_batches, config.error_file_path)

        start_time = time.time()
        combine_orphan_genes_stats(orphan_genes_internal_dir, orphan_genes_dir)
        step_post_processing_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} post-processing took {step_post_processing_time}.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


    # orthogroups_final
    step_number = f'{base_step_number}_{start_substep_number + 6}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_final'
    script_path = consts.SRC_DIR / 'steps' / 'add_orphans_to_orthogroups.py'
    final_orthogroups_dir_path, pipeline_step_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    final_orthogroups_file_path = final_orthogroups_dir_path / 'orthogroups.csv'
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Constructing final orthologs table...')

        if infer_orthogroups_config.add_orphan_genes_to_ogs:
            params = [orthogroups_file_path, final_orthogroups_file_path, f'--orphan_genes_dir {orphan_genes_internal_dir}']

            submit_mini_batch(logger, config, script_path, [params], pipeline_step_tmp_dir,
                              'add_orphans_to_orthogroups')
            wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, 1,
                             config.error_file_path)
        else:
            shutil.copy(orthogroups_file_path, final_orthogroups_file_path)
            logger.info(f'add_orphan_genes_to_ogs is False. Skipping adding orphans to orthogroups. Copied '
                        f'{orthogroups_file_path} to {final_orthogroups_file_path}.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return final_orthogroups_dir_path, orphan_genes_dir, start_substep_number + 6


def run_unified_mmseqs(logger, times_logger, config, infer_orthogroups_config, base_step_number):
    with open(infer_orthogroups_config.genomes_names_path) as f:
        strains_names = f.read().rstrip().split('\n')

    # 1.	mmseqs2_all_vs_all.py
    step_number = f'{base_step_number}_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_all_vs_all_analysis'
    script_path = consts.SRC_DIR / 'steps' / 'mmseqs_v2 '/ 'mmseqs2_all_vs_all.py'
    all_vs_all_output_dir, pipeline_step_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    m8_output_path = all_vs_all_output_dir / 'all_vs_all_reduced_columns.csv'
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        m8_raw_output_path = all_vs_all_output_dir / 'all_vs_all_raw.m8'

        cpus = min(consts.MMSEQS_BIG_DATASET_NUM_OF_CORES, infer_orthogroups_config.max_parallel_jobs)
        params = [infer_orthogroups_config.all_proteins_path,
                  all_vs_all_output_dir,
                  m8_raw_output_path,
                  f'--identity_cutoff {config.identity_cutoff / 100}',
                  f'--coverage_cutoff {config.coverage_cutoff / 100}',
                  f'--e_value_cutoff {config.e_value_cutoff}',
                  f'--sensitivity {config.sensitivity}',
                  f'--number_of_genomes {len(strains_names)}',
                  f'--cpus {cpus}']

        submit_mini_batch(logger, config, script_path, [params], pipeline_step_tmp_dir,
                          'mmseqs', num_of_cpus=cpus, memory=consts.MMSEQS_BIG_DATASET_REQUIRED_MEMORY_GB,
                          time_in_hours=consts.MMSEQS_JOB_TIME_LIMIT_HOURS)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, 1,
                         config.error_file_path)

        start_time = time.time()

        logger.info(f"Starting to read {m8_raw_output_path} and create a processed version of it...")
        m8_df = dd.read_csv(m8_raw_output_path, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER,
                            dtype=consts.MMSEQS_OUTPUT_COLUMNS_TYPES)
        m8_df = m8_df[m8_df['query'] != m8_df['target']]
        add_score_column_to_mmseqs_output(m8_df)
        m8_df['query_genome'] = m8_df['query'].str.split(':').str[0]
        m8_df['target_genome'] = m8_df['target'].str.split(':').str[0]
        m8_df = m8_df[['query', 'query_genome', 'target', 'target_genome', 'score']]
        m8_df.to_parquet(m8_output_path)  # Here I always use parquet since it's a huge file
        logger.info(f"{m8_output_path} was created successfully.")

        step_post_processing_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} post-processing took {step_post_processing_time}.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    m8_output_size_in_gb = get_directory_size_in_gb(m8_output_path)
    m8_output_parsing_memory = str(max(int(consts.DEFAULT_MEMORY_PER_JOB_GB), math.ceil(m8_output_size_in_gb * 10)))

    # 2.	extract_rbh_hits.py
    step_number = f'{base_step_number}_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_extract_rbh_hits'
    script_path = consts.SRC_DIR / 'steps' / 'mmseqs_v2' / 'extract_rbh_hits.py'
    orthologs_output_dir, pipeline_step_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    orthologs_scores_statistics_dir = orthologs_output_dir / 'scores_statistics'
    os.makedirs(orthologs_scores_statistics_dir, exist_ok=True)
    max_rbh_scores_parts_output_dir = orthologs_output_dir / 'max_rbh_scores_parts'
    os.makedirs(max_rbh_scores_parts_output_dir, exist_ok=True)
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        if len(strains_names) >= 2:
            rbh_inputs_dir = pipeline_step_tmp_dir / 'rbh_inputs_dir'
            os.makedirs(rbh_inputs_dir, exist_ok=True)
            job_index_to_pairs = defaultdict(list)
            for i, (genome1, genome2) in enumerate(itertools.combinations(strains_names, 2)):
                job_index = i % infer_orthogroups_config.max_parallel_jobs
                job_index_to_pairs[job_index].append((genome1, genome2))

            all_cmds_params = []
            for job_index, pairs in job_index_to_pairs.items():
                job_input_path = rbh_inputs_dir / f'job_{job_index}_pairs.txt'
                pairs_text = '\n'.join([f'{genome1}\t{genome2}' for genome1, genome2 in pairs])
                with open(job_input_path, 'w') as f:
                    f.write(pairs_text)

                params = [m8_output_path, job_input_path, orthologs_output_dir,
                          orthologs_scores_statistics_dir, max_rbh_scores_parts_output_dir]
                if config.use_parquet:
                    params.append('--use_parquet')
                if config.verbose:
                    params.append('--verbose')
                all_cmds_params.append(params)

            num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                          'rbh_analysis', memory=m8_output_parsing_memory,
                                          time_in_hours=consts.MMSEQS_JOB_TIME_LIMIT_HOURS)

            wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                             num_of_batches, config.error_file_path)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    if infer_orthogroups_config.skip_paralogs:
        return orthologs_output_dir, None, orthologs_scores_statistics_dir, None

    # 3. extract_paralogs.py
    step_number = f'{base_step_number}_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_paralogs'
    script_path = consts.SRC_DIR / 'steps' / 'mmseqs_v2' / 'extract_paralogs.py'
    paralogs_output_dir, pipeline_step_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    paralogs_scores_statistics_dir = paralogs_output_dir / 'scores_statistics'
    os.makedirs(paralogs_scores_statistics_dir, exist_ok=True)
    max_rbh_scores_unified_dir = paralogs_output_dir / 'max_rbh_scores_unified'
    os.makedirs(max_rbh_scores_unified_dir, exist_ok=True)
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Searching for paralogs in each genome')

        paralogs_inputs_dir = pipeline_step_tmp_dir / 'paralogs_inputs_dir'
        os.makedirs(paralogs_inputs_dir, exist_ok=True)
        job_index_to_genome = defaultdict(list)
        for i, genome in enumerate(strains_names):
            job_index = i % infer_orthogroups_config.max_parallel_jobs
            job_index_to_genome[job_index].append(genome)

        all_cmds_params = []
        for job_index, genomes in job_index_to_genome.items():
            job_input_path = paralogs_inputs_dir / f'job_{job_index}_genomes.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(genomes))

            single_cmd_params = [m8_output_path, job_input_path, max_rbh_scores_parts_output_dir,
                                 paralogs_output_dir, max_rbh_scores_unified_dir, paralogs_scores_statistics_dir]
            if config.use_parquet:
                single_cmd_params.append('--use_parquet')
            if config.verbose:
                single_cmd_params.append('--verbose')
            all_cmds_params.append(single_cmd_params)

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                      'paralogs_analysis', memory=m8_output_parsing_memory)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, config.error_file_path)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir


def run_non_unified_mmseqs_with_dbs(logger, times_logger, config, infer_orthogroups_config, base_step_number):
    with open(infer_orthogroups_config.genomes_names_path) as f:
        strains_names = f.read().rstrip().split('\n')

    # 1. mmseqs2_create_db
    step_number = f'{base_step_number}_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    script_path = consts.SRC_DIR / 'steps' / 'mmseqs_v4' / 'mmseqs2_create_db.py'
    step_name = f'{step_number}_mmseqs_dbs'
    mmseqs_dbs_output_dir, mmseqs_dbs_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Creating mmseqs dbs...')

        job_index_to_strain_names = defaultdict(list)

        for i, strain_name in enumerate(strains_names):
            job_index = i % infer_orthogroups_config.max_parallel_jobs
            job_index_to_strain_names[job_index].append(strain_name)

        jobs_inputs_dir = mmseqs_dbs_tmp_dir / 'job_inputs'
        os.makedirs(jobs_inputs_dir, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_strain_names in job_index_to_strain_names.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_strain_names))

            single_cmd_params = [job_input_path, infer_orthogroups_config.translated_orfs_dir, mmseqs_dbs_output_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, mmseqs_dbs_tmp_dir,
                                      'mmseqs_dbs')

        wait_for_results(logger, times_logger, step_name, mmseqs_dbs_tmp_dir,
                         num_of_batches, config.error_file_path)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 2.	mmseqs2_rbh.py
    step_number = f'{base_step_number}_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mmseqs_rbh'
    script_path = consts.SRC_DIR / 'steps' / 'mmseqs_v4' / 'mmseqs2_rbh.py'
    orthologs_output_dir, orthologs_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    orthologs_scores_statistics_dir = orthologs_output_dir / 'scores_statistics'
    os.makedirs(orthologs_scores_statistics_dir, exist_ok=True)
    max_rbh_scores_parts_output_dir = orthologs_output_dir / 'max_rbh_scores_parts'
    os.makedirs(max_rbh_scores_parts_output_dir, exist_ok=True)
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        if len(strains_names) >= 2:
            logger.info(f'Querying all VS all (using mmseqs2)...')
            temp_dir = orthologs_output_dir / 'temp'
            os.makedirs(temp_dir, exist_ok=True)

            job_index_to_strain_pairs = defaultdict(list)
            for i, (strain1_name, strain2_name) in enumerate(itertools.combinations(strains_names, 2)):
                job_index = i % infer_orthogroups_config.max_parallel_jobs
                job_index_to_strain_pairs[job_index].append((strain1_name, strain2_name))

            jobs_inputs_dir = orthologs_tmp_dir / 'job_inputs'
            os.makedirs(jobs_inputs_dir, exist_ok=True)

            all_cmds_params = []
            for job_index, job_strain_pairs in job_index_to_strain_pairs.items():
                job_input_path = jobs_inputs_dir / f'{job_index}.txt'
                with open(job_input_path, 'w') as f:
                    for strain1_name, strain2_name in job_strain_pairs:
                        f.write(f'{strain1_name}\t{strain2_name}\n')

                single_cmd_params = [job_input_path,
                                     mmseqs_dbs_output_dir,
                                     orthologs_output_dir,
                                     orthologs_scores_statistics_dir,
                                     max_rbh_scores_parts_output_dir,
                                     temp_dir,
                                     f'--identity_cutoff {config.identity_cutoff / 100}',
                                     f'--coverage_cutoff {config.coverage_cutoff / 100}',
                                     f'--e_value_cutoff {config.e_value_cutoff}',
                                     f'--sensitivity {config.sensitivity}',
                                     ]
                if config.use_parquet:
                    single_cmd_params.append('--use_parquet')
                all_cmds_params.append(single_cmd_params)

            num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, orthologs_tmp_dir,
                                          'rbh_analysis', time_in_hours=consts.MMSEQS_JOB_TIME_LIMIT_HOURS)

            wait_for_results(logger, times_logger, step_name, orthologs_tmp_dir,
                             num_of_batches, config.error_file_path)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    if infer_orthogroups_config.skip_paralogs:
        return orthologs_output_dir, None, orthologs_scores_statistics_dir, None

    # 3. mmseqs2_paralogs.py
    step_number = f'{base_step_number}_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mmseqs_paralogs'
    script_path = consts.SRC_DIR / 'steps' / 'mmseqs_v4' / 'mmseqs2_paralogs.py'
    paralogs_output_dir, paralogs_tmp_dir = prepare_directories(
        logger, infer_orthogroups_config.output_dir, infer_orthogroups_config.tmp_dir, step_name)
    paralogs_scores_statistics_dir = paralogs_output_dir / 'scores_statistics'
    os.makedirs(paralogs_scores_statistics_dir, exist_ok=True)
    max_rbh_scores_unified_dir = paralogs_output_dir / 'max_rbh_scores_unified'
    os.makedirs(max_rbh_scores_unified_dir, exist_ok=True)
    done_file_path = infer_orthogroups_config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Searching for paralogs in each genome')
        temp_dir = paralogs_output_dir / 'temp'
        os.makedirs(temp_dir, exist_ok=True)

        job_index_to_strain_names = defaultdict(list)

        for i, strain_name in enumerate(strains_names):
            job_index = i % infer_orthogroups_config.max_parallel_jobs
            job_index_to_strain_names[job_index].append(strain_name)

        jobs_inputs_dir = paralogs_tmp_dir / 'job_inputs'
        os.makedirs(jobs_inputs_dir, exist_ok=True)

        all_cmds_params = []

        for job_index, job_strains_names in job_index_to_strain_names.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_strains_names))

            single_cmd_params = [job_input_path,
                                 mmseqs_dbs_output_dir,
                                 max_rbh_scores_parts_output_dir,
                                 paralogs_output_dir,
                                 max_rbh_scores_unified_dir,
                                 paralogs_scores_statistics_dir,
                                 temp_dir,
                                 f'--identity_cutoff {config.identity_cutoff / 100}',
                                 f'--coverage_cutoff {config.coverage_cutoff / 100}',
                                 f'--e_value_cutoff {config.e_value_cutoff}',
                                 f'--sensitivity {config.sensitivity}',
                                 ]
            if config.use_parquet:
                single_cmd_params.append('--use_parquet')
            all_cmds_params.append(single_cmd_params)

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, paralogs_tmp_dir,
                                      'paralogs_analysis')

        wait_for_results(logger, times_logger, step_name, paralogs_tmp_dir, num_of_batches, config.error_file_path)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir
