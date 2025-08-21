import shutil
import sys
import math
from collections import defaultdict
import time
import traceback
import json
import subprocess
from datetime import timedelta

from Bio import SeqIO
import pandas as pd
from dataclasses import replace

from auxiliaries import consts
from auxiliaries.configuration import get_configuration
from auxiliaries.input_verifications import prepare_and_verify_input_data

from auxiliaries.general_utils import write_done_file, get_required_memory_gb_to_load_csv
from auxiliaries.main_utils import send_email_in_pipeline_end, submit_clean_folders_job, zip_results, \
    submit_clean_old_user_results_job, update_progressbar, define_intervals, add_results_to_final_dir, \
    calc_genomes_batch_size, initialize_progressbar, find_all_gap_sequences
from auxiliaries.run_step_utils import wait_for_results, prepare_directories, submit_job, submit_batch
from auxiliaries.logic_utils import (plot_genomes_histogram, combine_orphan_genes_stats,
                                     split_ogs_to_jobs_inputs_files_by_og_sizes, sort_orthogroups_df_and_rename_ogs)

from auxiliaries.infer_orthogroups_logic import infer_orthogroups
from flask.SharedConsts import State


def step_1_fix_input_files(logger, times_logger, config):
    # 01.	drop_plasmids_and_fix_frames.py
    step_number = '01'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_fix_input_files'
    script_path = consts.SRC_DIR / 'steps/drop_plasmids_and_fix_frames.py'
    filtered_inputs_dir, pipeline_step_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir,
                                                                     step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Filtering plasmids out...')

        job_index_to_fasta_files = defaultdict(list)

        for i, fasta_file in enumerate(config.data_path.iterdir()):
            job_index = i % config.max_parallel_jobs
            job_index_to_fasta_files[job_index].append(str(fasta_file))

        jobs_inputs_dir = pipeline_step_tmp_dir / 'job_inputs'
        jobs_inputs_dir.mkdir(parents=True, exist_ok=True)

        for job_index, job_fasta_files in job_index_to_fasta_files.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_fasta_files))

        script_params = [filtered_inputs_dir, f'--drop_plasmids {config.filter_out_plasmids}',
                         f'--fix_frames {config.inputs_fasta_type == "orfs"}']

        submit_batch(logger, config, script_path, script_params, jobs_inputs_dir, pipeline_step_tmp_dir,
                     'drop_plasmids')

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, config.error_file_path)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    if config.filter_out_plasmids:
        update_progressbar(config.progressbar_file_path, 'Filter out plasmids')

    return filtered_inputs_dir


def step_2_search_orfs(logger, times_logger, config):
    # 02_1.	search_orfs.py
    # Input: (1) an input path for a fasta file with contigs/full genome (2) an output file path
    #                (with a suffix as follows: i_genes.fasta. especially relevant for the wrapper).
    # Output: a fasta file where under each header, there’s a single gene.
    # Can be parallelized on cluster
    # Prodigal ignores newlines in the middle of a sequence, namely, >bac1\nAAA\nAA\nTTT >bac1\nAAAAATTT will be
    # analyzed identically.
    step_number = '02_1'  # using orfs_step_number for downstream analysis (extract promoter and gene sequences)
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orfs'
    script_path = consts.SRC_DIR / 'steps' / 'search_orfs.py'
    orfs_dir, orfs_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir, step_name)
    orfs_sequences_dir = orfs_dir / 'orfs_sequences'
    orfs_statistics_dir = orfs_dir / 'orfs_statistics'
    orfs_translated_dir = orfs_dir / 'orfs_translated'
    orfs_coordinates_dir = orfs_dir / 'orfs_coordinates'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Extracting ORFs...')

        orfs_sequences_dir.mkdir(parents=True, exist_ok=True)
        orfs_statistics_dir.mkdir(parents=True, exist_ok=True)
        orfs_translated_dir.mkdir(parents=True, exist_ok=True)
        orfs_coordinates_dir.mkdir(parents=True, exist_ok=True)

        job_index_to_fasta_files = defaultdict(list)

        for i, fasta_file in enumerate(config.data_path.iterdir()):
            job_index = i % config.max_parallel_jobs
            job_index_to_fasta_files[job_index].append(str(fasta_file))

        jobs_inputs_dir = orfs_tmp_dir / 'job_inputs'
        jobs_inputs_dir.mkdir(parents=True, exist_ok=True)

        for job_index, job_fasta_files in job_index_to_fasta_files.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_fasta_files))

        script_params = [orfs_sequences_dir, orfs_statistics_dir, orfs_translated_dir, orfs_coordinates_dir, config.inputs_fasta_type]
        submit_batch(logger, config, script_path, script_params, jobs_inputs_dir, orfs_tmp_dir, 'search_orfs')

        wait_for_results(logger, times_logger, step_name, orfs_tmp_dir, config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orfs_sequences_dir, config.final_output_dir)
            add_results_to_final_dir(logger, orfs_translated_dir, config.final_output_dir)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 02_2.	plot_orfs_statistics
    step_number = '02_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orfs_plots'
    orfs_plots_path, pipeline_step_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir,
                                                                 step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Concatenating orfs counts and gc contents...')
        start_time = time.time()

        gc_content = {}
        orf_count = {}
        for orfs_file in orfs_statistics_dir.iterdir():
            strain_name = orfs_file.stem
            with open(orfs_file, 'r') as fp:
                orfs_statistics = json.load(fp)

            gc_content[strain_name] = orfs_statistics['gc_content']
            orf_count[strain_name] = orfs_statistics['orfs_count']

        plot_genomes_histogram(orf_count, orfs_plots_path, 'orfs_counts', 'ORFs count', 'ORFs Count per genome')
        plot_genomes_histogram(gc_content, orfs_plots_path, 'orfs_gc_content', 'GC Content (of ORFs)',
                               'GC Content per genome')

        step_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} took {step_time}.')

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orfs_plots_path, config.final_output_dir)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 02_3.  Aggregate all protein fasta files to one file
    step_number = '02_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_concat_orfs'
    concat_orfs_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir,
                                                                      step_name)
    all_orfs_fasta_path = concat_orfs_dir_path / 'all_orfs.fna'
    all_proteins_fasta_path = concat_orfs_dir_path / 'all_proteomes.faa'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Concatenating ORFs sequences...')
        start_time = time.time()

        cmd = f"cat {orfs_sequences_dir / '*'} > {all_orfs_fasta_path}"
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

        cmd = f"cat {orfs_translated_dir / '*'} > {all_proteins_fasta_path}"
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

        step_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} took {step_time}.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    update_progressbar(config.progressbar_file_path, 'Predict and translate ORFs')
    update_progressbar(config.progressbar_file_path, 'Translate ORFs')

    return orfs_sequences_dir, orfs_translated_dir, all_orfs_fasta_path, all_proteins_fasta_path, orfs_coordinates_dir


def step_3_analyze_genome_completeness(logger, times_logger, config, translated_orfs_dir_path):
    # 03.  assessing_genomes_completeness.py
    # Input: path to faa file
    # Output: calculate proteome completeness based on a dataset of profile-HMMs that represent core bacterial genes.
    step_number = '03'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_genomes_completeness'
    script_path = consts.SRC_DIR / 'steps' / 'assessing_genome_completeness.py'
    genome_completeness_dir_path, genome_completeness_tmp_dir = prepare_directories(logger, config.steps_results_dir,
                                                                                    config.tmp_dir, step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Calculating genomes completeness...')

        job_index_to_fasta_files = defaultdict(list)

        for i, fasta_file in enumerate(translated_orfs_dir_path.iterdir()):
            job_index = i % config.max_parallel_jobs
            job_index_to_fasta_files[job_index].append(str(fasta_file))

        jobs_inputs_dir = genome_completeness_tmp_dir / 'job_inputs'
        jobs_inputs_dir.mkdir(parents=True, exist_ok=True)
        genomes_output_dir_path = genome_completeness_tmp_dir / 'individual_proteomes_outputs'
        genomes_output_dir_path.mkdir(parents=True, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_fasta_files in job_index_to_fasta_files.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_fasta_files))

            single_cmd_params = [job_input_path, genomes_output_dir_path]
            all_cmds_params.append(single_cmd_params)

        submit_batch(logger, config, script_path, all_cmds_params, genome_completeness_tmp_dir,
                                      'genomes_completeness')

        wait_for_results(logger, times_logger, step_name, genome_completeness_tmp_dir, config.error_file_path)

        # Aggregate results of step
        start_time = time.time()

        genomes_completeness_scores = {}
        for genome_dir in genomes_output_dir_path.iterdir():
            genome_score_path = genome_dir / 'result.txt'
            with open(genome_score_path, 'r') as fp:
                genomes_completeness_scores[genome_dir.name] = float(fp.read().strip())

        plot_genomes_histogram(genomes_completeness_scores, genome_completeness_dir_path, 'genomes_completeness',
                               'Genome BUSCO completeness score', 'Genome BUSCO completeness score')

        shutil.rmtree(genomes_output_dir_path, ignore_errors=True)

        aggregation_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} post-processing took {aggregation_time}.')

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, genome_completeness_dir_path, config.final_output_dir)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    update_progressbar(config.progressbar_file_path, 'Calculate genomes completeness')


def step_5_orthogroups_inference(logger, times_logger, config, genomes_names, translated_orfs_dir, all_proteins_path,
                                 orfs_coordinates_dir):
    genomes_batch_size = calc_genomes_batch_size(logger, config, len(genomes_names))
    if len(genomes_names) < config.min_num_of_genomes_to_optimize_orthogroups_inference or \
            len(genomes_names) < genomes_batch_size * 2 or config.always_run_full_orthogroups_inference:
        final_orthogroups_dir_path, orphan_genes_dir, final_substep_number = infer_orthogroups(
            logger, times_logger, config, '05', translated_orfs_dir, all_proteins_path)
        final_step_number = '05'
    else:
        final_orthogroups_dir_path, orphan_genes_dir, final_substep_number = step_5_6_approximate_orthogroups_inference(
            logger, times_logger, config, translated_orfs_dir, genomes_batch_size)
        final_step_number = '06'

    # Sort orthogroups table by coordinates
    step_number = f'{final_step_number}_{final_substep_number + 1}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_sort_orthogroups_by_coordinates'
    sorted_orthogroups_dir, sorted_orthogroups_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, step_name)
    sorted_orthogroups_file_path = sorted_orthogroups_dir / 'orthogroups.csv'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        start_time = time.time()

        sort_orthogroups_df_and_rename_ogs(logger, final_orthogroups_dir_path / 'orthogroups.csv',
                                           orfs_coordinates_dir, sorted_orthogroups_file_path)

        step_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} took {step_time}.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    if not config.do_not_copy_outputs_to_final_results_dir:
        add_results_to_final_dir(logger, orphan_genes_dir, config.final_output_dir)

    if not config.do_not_copy_outputs_to_final_results_dir:
        add_results_to_final_dir(logger, sorted_orthogroups_dir, config.final_output_dir)

    update_progressbar(config.progressbar_file_path, 'Infer orthogroups')
    update_progressbar(config.progressbar_file_path, 'Find orphan genes')

    return sorted_orthogroups_file_path


def step_5_6_approximate_orthogroups_inference(logger, times_logger, config, translated_orfs_dir, genomes_batch_size):
    with open(config.genomes_names_path, 'r') as genomes_names_fp:
        genomes_names = genomes_names_fp.read().split('\n')

    # 05.	subsets_inference
    step_number = '05'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_subsets_inference'
    script_path = consts.SRC_DIR / 'steps' / 'infer_orthogroups_on_genomes_batch.py'
    inference_dir_path, inference_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir,
                                                                step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    done_dir_path = config.done_files_dir / step_name
    if not done_file_path.exists():
        genomes_batches = define_intervals(len(genomes_names), genomes_batch_size)
        logger.info(f'Infer orthogroups on {len(genomes_names)} genomes in {len(genomes_batches)} batches, '
                    f'using batch size of {genomes_batch_size}')

        job_index_to_batches = defaultdict(list)

        for batch_id, (start_index, end_index_inclusive) in enumerate(genomes_batches):
            subset_genome_names = genomes_names[start_index:end_index_inclusive + 1]

            subset_output_dir = inference_dir_path / f'subset_{batch_id}'
            subset_tmp_dir = inference_tmp_dir / f'subset_{batch_id}'
            subset_done_dir = done_dir_path / f'subset_{batch_id}'
            subset_output_dir.mkdir(parents=True, exist_ok=True)
            subset_tmp_dir.mkdir(parents=True, exist_ok=True)
            subset_done_dir.mkdir(parents=True, exist_ok=True)

            subset_genomes_names_path = subset_output_dir / f'genomes_names.txt'
            with open(subset_genomes_names_path, 'w') as subset_genomes_names_fp:
                subset_genomes_names_fp.write('\n'.join(subset_genome_names))

            config_for_subset_inference = replace(
                config, genomes_names_path=subset_genomes_names_path, steps_results_dir=subset_output_dir,
                tmp_dir=subset_tmp_dir, done_files_dir=subset_done_dir,
                max_parallel_jobs=max(1, config.max_parallel_jobs // len(genomes_batches)),
                add_orphan_genes_to_ogs=True)
            subset_config_path = subset_tmp_dir / 'config.csv'
            config_for_subset_inference.to_csv(subset_config_path)

            job_index = batch_id % config.max_parallel_jobs
            job_index_to_batches[job_index].append(f'{batch_id}\t{subset_config_path}')

        jobs_inputs_dir = inference_tmp_dir / 'job_inputs'
        jobs_inputs_dir.mkdir(parents=True, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_batches in job_index_to_batches.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_batches))

            params = [job_input_path, step_number, translated_orfs_dir]
            all_cmds_params.append(params)

        submit_batch(logger, config, script_path, all_cmds_params, inference_tmp_dir,
                                      'infer_orthogroups', time_in_hours=config.infer_orthogroups_time_limit)

        wait_for_results(logger, times_logger, step_name, inference_tmp_dir, config.error_file_path, recursive_step=True)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 06_0.	aggregate_sub_orthogroups
    step_number = '06_0'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_aggregate_sub_orthogroups'
    aggregate_orthogroups_dir_path, aggregate_orthogroups_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, step_name)
    pseudo_genomes_strains_names_path = aggregate_orthogroups_dir_path / 'pseudo_genomes_names.txt'
    pseudo_genomes_dir_path = aggregate_orthogroups_dir_path / 'pseudo_genomes'
    all_pseudo_genomes_path = aggregate_orthogroups_dir_path / f'all_pseudo_proteomes.faa'
    sub_orthogroups_dir_path = aggregate_orthogroups_dir_path / 'sub_orthogroups'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Collection psuedo genomes and sub-orthogroups...')
        start_time = time.time()

        pseudo_genomes_dir_path.mkdir(parents=True, exist_ok=True)
        sub_orthogroups_dir_path.mkdir(parents=True, exist_ok=True)
        pseudo_genomes_names = []

        for batch_dir_path in inference_dir_path.iterdir():
            if not batch_dir_path.is_dir():
                continue
            batch_id = batch_dir_path.name.split('_')[1]

            orthogroups_with_representative_path = batch_dir_path / '05_12_pseudo_genome' / f'orthogroups_with_representative_{batch_id}.csv'
            pseudo_genome_file_path = batch_dir_path / '05_12_pseudo_genome' / f'pseudo_genome_{batch_id}.faa'

            shutil.copy(orthogroups_with_representative_path, sub_orthogroups_dir_path)
            shutil.copy(pseudo_genome_file_path, pseudo_genomes_dir_path)
            pseudo_genomes_names.append(f'pseudo_genome_{batch_id}')

        cmd = f"cat {pseudo_genomes_dir_path / '*'} > {all_pseudo_genomes_path}"
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

        with open(pseudo_genomes_strains_names_path, 'w') as pseudo_genomes_strains_names_fp:
            pseudo_genomes_strains_names_fp.write('\n'.join(pseudo_genomes_names))

        step_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} took {step_time}.')
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 06. infer pseudo orthogroups
    config_for_pseudo_orthogroups_inference = replace(
        config, genomes_names_path=pseudo_genomes_strains_names_path, add_orphan_genes_to_ogs=True)
    pseudo_orthogroups_dir_path, _, final_substep_number = infer_orthogroups(
        logger, times_logger, config_for_pseudo_orthogroups_inference, '06', pseudo_genomes_dir_path,
        all_pseudo_genomes_path, skip_paralogs=True)

    pseudo_orthogroups_file_path = pseudo_orthogroups_dir_path / 'orthogroups.csv'

    # 06_11. merged_suborthogroups
    step_number = f'06_{final_substep_number + 1}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_merged_suborthogroups'
    script_path = consts.SRC_DIR / 'steps' / 'merge_sub_orthogroups.py'
    merged_orthogroups_dir_path, merged_orthogroups_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, step_name)
    merged_orthogroups_file_path = merged_orthogroups_dir_path / 'orthogroups.csv'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Merge sub-orthogroups...')

        params = [pseudo_orthogroups_file_path, sub_orthogroups_dir_path, merged_orthogroups_file_path]

        submit_job(logger, config, script_path, params, merged_orthogroups_tmp_dir,
                          'merge_sub_orthogroups', memory=config.merge_sub_orthogroups_memory_gb)
        wait_for_results(logger, times_logger, step_name, merged_orthogroups_tmp_dir, config.error_file_path)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 06.12 extract_orphan_genes_from_full_orthogroups.py
    step_number = f'06_{final_substep_number + 2}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orphan_genes_from_orthogroups'
    script_path = consts.SRC_DIR / 'steps' / 'extract_orphan_genes_from_full_orthogroups.py'
    orphan_genes_dir, orphans_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir, step_name)
    orphan_genes_internal_dir = orphan_genes_dir / 'orphans_lists_per_genome'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Extracting orphan genes...')
        orphan_genes_internal_dir.mkdir(parents=True, exist_ok=True)

        job_index_to_genome_names = defaultdict(list)

        for i, genome_name in enumerate(genomes_names):
            job_index = i % config.max_parallel_jobs
            job_index_to_genome_names[job_index].append(genome_name)

        jobs_inputs_dir = orphans_tmp_dir / 'job_inputs'
        jobs_inputs_dir.mkdir(parents=True, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_genome_names in job_index_to_genome_names.items():
            job_input_path = jobs_inputs_dir / f'{job_index}.txt'
            with open(job_input_path, 'w') as f:
                f.write('\n'.join(job_genome_names))

            single_cmd_params = [job_input_path, merged_orthogroups_file_path, orphan_genes_internal_dir]
            all_cmds_params.append(single_cmd_params)

        memory_gb = max(config.job_default_memory_gb, get_required_memory_gb_to_load_csv(merged_orthogroups_file_path))
        submit_batch(logger, config, script_path, all_cmds_params, orphans_tmp_dir,
                                      'extract_orphans_from_orthogroups', memory=memory_gb)

        wait_for_results(logger, times_logger, step_name, orphans_tmp_dir, config.error_file_path)

        start_time = time.time()
        combine_orphan_genes_stats(orphan_genes_internal_dir, orphan_genes_dir)
        step_post_processing_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} post-processing took {step_post_processing_time}.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 06.13 final_orthogroups_table
    step_number = f'06_{final_substep_number + 3}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_final'
    final_orthogroups_dir, final_orthogroups_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, step_name)
    final_orthogroups_file_path = final_orthogroups_dir / 'orthogroups.csv'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        start_time = time.time()

        if config.add_orphan_genes_to_ogs:
            shutil.copy(merged_orthogroups_file_path, final_orthogroups_file_path)
            logger.info(f'add_orphan_genes_to_ogs is True. Copied {merged_orthogroups_file_path} to '
                        f'{final_orthogroups_file_path} since it already contains orphans.')
        else:
            orthogroups_df = pd.read_csv(merged_orthogroups_file_path, index_col='OG_name', dtype=str)
            orthogroups_df = orthogroups_df[~(
                    (orthogroups_df.count(axis=1) == 1) &
                    ~(orthogroups_df.apply(lambda row: row.dropna().iloc[0].__contains__(';'), axis=1))
            )]
            orthogroups_df.reset_index(drop=True, inplace=True)
            orthogroups_df['OG_name'] = [f'OG_{i}' for i in range(len(orthogroups_df.index))]
            orthogroups_df.set_index('OG_name', inplace=True)
            orthogroups_df.to_csv(final_orthogroups_file_path)
            logger.info(
                f'add_orphan_genes_to_ogs is False. Removed single orphan genes from {merged_orthogroups_file_path} '
                f'and saved to {final_orthogroups_file_path}')

        step_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} took {step_time}.')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return final_orthogroups_dir, orphan_genes_dir, final_substep_number + 3


def step_7_orthologs_table_variations(logger, times_logger, config, final_orthogroups_file_path, orfs_coordinates_dir):
    # 07_1.	create_orthoxml.py
    step_number = '07_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthoxml'
    script_path = consts.SRC_DIR / 'steps' / 'create_orthoxml.py'
    orthoxml_dir_path, orthoxml_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir, step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Creating orthoxml...')

        params = [final_orthogroups_file_path, orthoxml_dir_path, f'--qfo_benchmark {config.qfo_benchmark}']

        submit_job(logger, config, script_path, params, orthoxml_tmp_dir,
                          'orthoxml', memory=config.orthoxml_memory_gb)
        wait_for_results(logger, times_logger, step_name, orthoxml_tmp_dir, config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orthoxml_dir_path, config.final_output_dir)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 07_2.	orthogroups_visualizations.py
    step_number = '07_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_visualizations'
    script_path = consts.SRC_DIR / 'steps' / 'orthogroups_visualizations.py'
    visualizations_dir_path, visualizations_tmp_dir = prepare_directories(logger, config.steps_results_dir,
                                                                          config.tmp_dir, step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Creating orthogroups visualizations...')

        params = [final_orthogroups_file_path, orfs_coordinates_dir, visualizations_dir_path]

        submit_job(logger, config, script_path, params, visualizations_tmp_dir,
                          'orthogroups_visualizations', memory=config.orthogroups_visualizations_memory_gb)
        wait_for_results(logger, times_logger, step_name, visualizations_tmp_dir, config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, visualizations_dir_path, config.final_output_dir)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return visualizations_dir_path


def step_8_build_orthologous_groups_fastas(logger, times_logger, config, all_orfs_fasta_path,
                                           all_proteins_fasta_path, final_orthologs_table_file_path):
    # 08_1.	extract_orfs.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    step_number = '08_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_fasta'
    script_path = consts.SRC_DIR / 'steps' / 'extract_orfs.py'
    orthogroups_fasta_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, config.steps_results_dir,
                                                                            config.tmp_dir, step_name)
    orthogroups_dna_dir_path = orthogroups_fasta_dir_path / 'orthogroups_dna'
    orthologs_aa_dir_path = orthogroups_fasta_dir_path / 'orthogroups_aa'
    orthogroups_aa_msa_dir_path = orthogroups_fasta_dir_path / 'orthogroups_aa_msa'
    orthogroups_induced_dna_msa_dir_path = orthogroups_fasta_dir_path / 'orthogroups_induced_dna_msa'
    orthogroups_aa_consensus_dir_path = orthogroups_fasta_dir_path / 'orthogroups_aa_consensus'

    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')
        start_time = time.time()

        orthogroups_dna_dir_path.mkdir(parents=True, exist_ok=True)
        orthologs_aa_dir_path.mkdir(parents=True, exist_ok=True)
        orthogroups_aa_msa_dir_path.mkdir(parents=True, exist_ok=True)
        orthogroups_induced_dna_msa_dir_path.mkdir(parents=True, exist_ok=True)
        orthogroups_aa_consensus_dir_path.mkdir(parents=True, exist_ok=True)

        orthogroups_df = pd.read_csv(final_orthologs_table_file_path, dtype=str)
        job_paths = split_ogs_to_jobs_inputs_files_by_og_sizes(orthogroups_df, pipeline_step_tmp_dir,
                                                               config.max_parallel_jobs)
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_path in job_paths:
            single_cmd_params = [all_orfs_fasta_path,
                                 all_proteins_fasta_path,
                                 final_orthologs_table_file_path,
                                 job_path,
                                 orthogroups_dna_dir_path,
                                 orthologs_aa_dir_path,
                                 orthogroups_aa_msa_dir_path,
                                 orthogroups_induced_dna_msa_dir_path,
                                 orthogroups_aa_consensus_dir_path if config.kegg_optimization_mode == 'consensus_of_og' else None]
            all_cmds_params.append(single_cmd_params)

        orfs_size_gb = (all_orfs_fasta_path.stat().st_size + all_proteins_fasta_path.stat().st_size) / 1024 ** 3
        memory_gb = max(config.job_default_memory_gb,
                        math.ceil(orfs_size_gb * 2) + get_required_memory_gb_to_load_csv(final_orthologs_table_file_path))

        step_pre_processing_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} pre-processing time took {step_pre_processing_time}.')

        submit_batch(logger, config, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                      'orfs_extraction', memory=memory_gb)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orthogroups_dna_dir_path, config.final_output_dir)
            add_results_to_final_dir(logger, orthologs_aa_dir_path, config.final_output_dir)
            add_results_to_final_dir(logger, orthogroups_aa_msa_dir_path, config.final_output_dir)
            add_results_to_final_dir(logger, orthogroups_induced_dna_msa_dir_path, config.final_output_dir)

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    update_progressbar(config.progressbar_file_path, 'Prepare orthogroups fasta files')

    return orthogroups_dna_dir_path, orthologs_aa_dir_path, orthogroups_aa_msa_dir_path, \
        orthogroups_induced_dna_msa_dir_path, orthogroups_aa_consensus_dir_path


def step_9_extract_core_genome_and_core_proteome(logger, times_logger, config, aa_alignments_path, dna_alignments_path):
    # 09_1.	extract aligned_core_proteome.py
    step_number = '09_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    core_proteome_step_name = f'{step_number}_aligned_core_proteome'
    script_path = consts.SRC_DIR / 'steps' / 'extract_core_genome.py'
    aligned_core_proteome_path, aligned_core_proteome_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, core_proteome_step_name)
    aligned_core_proteome_file_path = aligned_core_proteome_path / 'aligned_core_proteome.fas'
    core_proteome_length_file_path = aligned_core_proteome_path / 'core_length.txt'
    core_proteome_done_file_path = config.done_files_dir / f'{core_proteome_step_name}.txt'
    if not core_proteome_done_file_path.exists():
        logger.info('Extracting aligned core proteome...')

        params = [aa_alignments_path, config.genomes_names_path,
                  aligned_core_proteome_file_path,
                  core_proteome_length_file_path,
                  f'--core_minimal_percentage {config.core_minimal_percentage}']  # how many members induce a core group?
        submit_job(logger, config, script_path, params, aligned_core_proteome_tmp_dir,
                          'core_proteome')

    else:
        logger.info(f'done file {core_proteome_done_file_path} already exists. Skipping step...')

    # 09_2.      extract_aligned_core_genome
    step_number = '09_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    core_genome_step_name = f'{step_number}_aligned_core_genome'
    script_path = consts.SRC_DIR / 'steps' / 'extract_core_genome.py'
    aligned_core_genome_path, aligned_core_genome_tmp_dir = prepare_directories(logger, config.steps_results_dir,
                                                                                config.tmp_dir, core_genome_step_name)
    core_genome_done_file_path = config.done_files_dir / f'{core_genome_step_name}.txt'
    if not core_genome_done_file_path.exists():
        logger.info('Extracting aligned core genome...')

        aligned_core_genome_file_path = aligned_core_genome_path / 'aligned_core_genome.fas'
        core_genome_length_file_path = aligned_core_genome_path / 'core_length.txt'

        params = [dna_alignments_path, config.genomes_names_path,
                  aligned_core_genome_file_path,
                  core_genome_length_file_path,
                  f'--core_minimal_percentage {config.core_minimal_percentage}']  # how many members induce a core group?
        submit_job(logger, config, script_path, params, aligned_core_genome_tmp_dir, 'core_genome')

    else:
        logger.info(f'done file {core_genome_done_file_path} already exists. Skipping step...')

    # 09_3.	extract aligned_core_proteome.py (for phylogeny reconstruction)
    step_number = '09_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    core_proteome_reduced_step_name = f'{step_number}_aligned_core_proteome_reduced'
    script_path = consts.SRC_DIR / 'steps' / 'extract_core_genome.py'
    aligned_core_proteome_reduced_path, aligned_core_proteome_reduced_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, core_proteome_reduced_step_name)
    aligned_core_proteome_reduced_file_path = aligned_core_proteome_reduced_path / 'aligned_core_proteome.fas'
    core_proteome_reduced_length_file_path = aligned_core_proteome_reduced_path / 'core_length.txt'
    core_proteome_reduced_done_file_path = config.done_files_dir / f'{core_proteome_reduced_step_name}.txt'
    if not core_proteome_reduced_done_file_path.exists():
        if config.max_number_of_core_ogs_for_phylogeny != -1:
            logger.info('Extracting aligned core proteome for phylogeny reconstruction...')

            params = [aa_alignments_path, config.genomes_names_path,
                      aligned_core_proteome_reduced_file_path,
                      core_proteome_reduced_length_file_path,
                      f'--core_minimal_percentage {config.core_minimal_percentage}',  # how many members induce a core group?
                      f'--max_number_of_ogs {config.max_number_of_core_ogs_for_phylogeny}']
            submit_job(logger, config, script_path, params, aligned_core_proteome_reduced_tmp_dir,
                              'core_proteome')
        else:
            logger.info('max_number_of_core_ogs_for_phylogeny is -1. Skipping step...')
            write_done_file(logger, core_proteome_reduced_done_file_path)
    else:
        logger.info(f'done file {core_proteome_reduced_done_file_path} already exists. Skipping step...')

    # Wait for all results here (since the steps aren't dependent on each other)
    if not core_proteome_done_file_path.exists():
        wait_for_results(logger, times_logger, core_proteome_step_name, aligned_core_proteome_tmp_dir,
                         config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, aligned_core_proteome_path, config.final_output_dir)
        write_done_file(logger, core_proteome_done_file_path)

    if not core_genome_done_file_path.exists():
        wait_for_results(logger, times_logger, core_genome_step_name, aligned_core_genome_tmp_dir,
                         config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, aligned_core_genome_path, config.final_output_dir)
        write_done_file(logger, core_genome_done_file_path)

    if not core_proteome_reduced_done_file_path.exists():
        wait_for_results(logger, times_logger, core_proteome_reduced_step_name, aligned_core_proteome_reduced_tmp_dir,
                         config.error_file_path)
        write_done_file(logger, core_proteome_reduced_done_file_path)

    # if step 09_3 was skipped, meaning we don't limit the number of core OGs for phylogeny reconstruction
    if config.max_number_of_core_ogs_for_phylogeny == -1:
        aligned_core_proteome_reduced_file_path = aligned_core_proteome_file_path
        core_proteome_reduced_length_file_path = core_proteome_length_file_path

    with open(core_proteome_reduced_length_file_path, 'r') as fp:
        core_proteome_reduced_length = int(fp.read().strip())

    update_progressbar(config.progressbar_file_path, 'Infer core genome and proteome')
    return aligned_core_proteome_reduced_file_path, core_proteome_reduced_length


def step_11_phylogeny(logger, times_logger, config, aligned_core_proteome_file_path, genomes_names,
                      core_proteome_length):
    # 11_1.   ani.py
    step_number = '11_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    ani_step_name = f'{step_number}_ani'
    script_path = consts.SRC_DIR / 'steps' / 'ani.py'
    ani_output_dir, ani_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir, ani_step_name)
    ani_done_file_path = config.done_files_dir / f'{ani_step_name}.txt'
    if not ani_done_file_path.exists():
        logger.info('Calculating ANI values...')
        ani_tmp_files = ani_tmp_dir / 'temp_results'
        ani_tmp_files.mkdir(parents=True, exist_ok=True)

        genomes_paths = [str(genome_path) for genome_path in config.data_path.iterdir()]
        genomes_list_path = ani_tmp_files / 'genomes_list.txt'
        with open(genomes_list_path, 'w') as genomes_list_file:
            genomes_list_file.write('\n'.join(genomes_paths))

        single_cmd_params = [genomes_list_path, ani_output_dir]

        submit_job(logger, config, script_path, single_cmd_params, ani_tmp_dir,
                          'ANI', num_of_cpus=config.ani_cpus, memory=config.ani_memory_gb)

    else:
        logger.info(f'done file {ani_done_file_path} already exists. Skipping step...')

    # 11_2.	reconstruct_species_phylogeny.py
    step_number = '11_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    phylogeny_step_name = f'{step_number}_species_phylogeny'
    script_path = consts.SRC_DIR / 'steps' / 'reconstruct_species_phylogeny.py'
    phylogeny_path, phylogeny_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir,
                                                            phylogeny_step_name)
    phylogeny_done_file_path = config.done_files_dir / f'{phylogeny_step_name}.txt'
    if not phylogeny_done_file_path.exists():
        output_tree_path = phylogeny_path / 'final_species_tree.newick'

        all_gap_seqs = find_all_gap_sequences(aligned_core_proteome_file_path)

        if len(genomes_names) < 3:
            error_text = 'Not enough genomes for phylogeny reconstruction. At least 3 genomes are required.'
        elif core_proteome_length == 0:
            error_text = ('Core proteome is empty. No core genes were found in the input genomes, hence skipping '
                          'phylogeny reconstruction.')
        elif all_gap_seqs:
            error_text = (f'Skipping phylogeny reconstruction since the genomes {",".join(sorted(all_gap_seqs))} have '
                          f'no genes at all in the core proteome.')
        else:
            error_text = None

        if error_text:
            logger.info(error_text)
            with open(output_tree_path, 'w') as output_tree_fp:
                output_tree_fp.write(error_text)
            output_tree_image_path = output_tree_path.with_suffix('.png')
            open(output_tree_image_path, 'w').close()
        else:
            logger.info('Reconstructing species phylogeny...')

            params = [aligned_core_proteome_file_path,
                      output_tree_path,
                      phylogeny_tmp_dir,
                      f'--bootstrap {config.bootstrap}']
            if config.outgroup:
                if config.outgroup in genomes_names:
                    params += [f'--outgroup {config.outgroup}']
                else:
                    logger.info(f'Outgroup {config.outgroup} was specified but it is not one of the input species:\n'
                                f'{",".join(sorted(genomes_names))}\nAn unrooted tree is going to be reconstructed')

            # Needed to avoid an error in drawing the tree (since the default xdg_runtime_dir is sometimes not writable)
            xdg_runtime_dir = phylogeny_tmp_dir / 'xdg_runtime_dir'
            xdg_runtime_dir.mkdir(parents=True, exist_ok=True)
            xdg_runtime_dir.chmod(0o700)

            submit_job(logger, config, script_path, params, phylogeny_tmp_dir,
                              'tree_reconstruction', num_of_cpus=config.phylogeny_cpus,
                       memory=config.phylogeny_memory_gb,
                       # Needed to avoid an error in drawing the tree. Taken from: https://github.com/NVlabs/instant-ngp/discussions/300
                       environment_variables_to_change_before_script={'QT_QPA_PLATFORM': 'offscreen', 'XDG_RUNTIME_DIR': xdg_runtime_dir},
                       time_in_hours=config.phylogeny_time_limit)

            # wait for the phylogenetic tree here
            wait_for_results(logger, times_logger, phylogeny_step_name, phylogeny_tmp_dir, config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, phylogeny_path, config.final_output_dir)
        write_done_file(logger, phylogeny_done_file_path)
    else:
        logger.info(f'done file {phylogeny_done_file_path} already exists. Skipping step...')

    # Wait for ANI here, such that ANI and phylogeny will run in parallel
    if not ani_done_file_path.exists():
        wait_for_results(logger, times_logger, ani_step_name, ani_tmp_dir, config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, ani_output_dir, config.final_output_dir)
        write_done_file(logger, ani_done_file_path)

    update_progressbar(config.progressbar_file_path, 'Reconstruct species phylogeny')
    update_progressbar(config.progressbar_file_path, 'Calculate ANI (Average Nucleotide Identity)')


def step_12_orthogroups_annotations(logger, times_logger, config, orfs_dir,
                                    orthologs_dna_sequences_dir_path, ogs_aa_sequences_dir_path,
                                    ogs_aa_consensus_dir_path, final_orthologs_table_file_path,
                                    orthogroups_visualizations_dir_path):
    # 12_1.  codon_bias.py
    # Input: ORF dir and OG dir
    # Output: W_vector for each genome, CAI for each OG
    step_number = '12_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    codon_bias_step_name = f'{step_number}_codon_bias'
    script_path = consts.SRC_DIR / 'steps' / 'codon_bias.py'
    codon_bias_output_dir_path, codon_bias_tmp_dir = prepare_directories(logger, config.steps_results_dir,
                                                                         config.tmp_dir, codon_bias_step_name)
    cai_table_path = codon_bias_output_dir_path / 'CAI_table.csv'
    codon_bias_done_file_path = config.done_files_dir / f'{codon_bias_step_name}.txt'
    if not codon_bias_done_file_path.exists():
        logger.info('Analyzing codon bias...')

        params = [
            orfs_dir,
            orthologs_dna_sequences_dir_path,
            codon_bias_output_dir_path,
            cai_table_path
        ]
        submit_job(logger, config, script_path, params, codon_bias_tmp_dir,
                          'codon_bias', num_of_cpus=config.codon_bias_cpus, memory=config.codon_bias_memory)

    else:
        logger.info(f'done file {codon_bias_done_file_path} already exists. Skipping step...')

    # 12_2.  kegg_annotation.py
    # Input: OG aa dir
    step_number = '12_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    kegg_step_name = f'{step_number}_kegg'
    script_path = consts.SRC_DIR / 'steps' / 'kegg_annotation.py'
    kegg_output_dir_path, kegg_tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir,
                                                             kegg_step_name)
    kegg_table_path = kegg_output_dir_path / 'og_kegg.csv'
    kegg_done_file_path = config.done_files_dir / f'{kegg_step_name}.txt'
    if not kegg_done_file_path.exists():
        logger.info('Annotation with KEGG Orthology...')

        params = [
            ogs_aa_sequences_dir_path if config.kegg_optimization_mode != 'consensus_of_og' else None,
            ogs_aa_consensus_dir_path if config.kegg_optimization_mode == 'consensus_of_og' else None,
            final_orthologs_table_file_path,
            kegg_output_dir_path,
            kegg_table_path,
            f'--optimization_mode {config.kegg_optimization_mode}'
        ]
        submit_job(logger, config, script_path, params, kegg_tmp_dir,
                          'kegg', num_of_cpus=config.kegg_cpus, memory=config.kegg_memory_gb)

    else:
        logger.info(f'done file {kegg_done_file_path} already exists. Skipping step...')

    # Wait for results here (since the steps aren't dependent on each other)
    if not codon_bias_done_file_path.exists():
        wait_for_results(logger, times_logger, codon_bias_step_name, codon_bias_tmp_dir, config.error_file_path)

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, codon_bias_output_dir_path, config.final_output_dir)
        write_done_file(logger, codon_bias_done_file_path)

    if not kegg_done_file_path.exists():
        wait_for_results(logger, times_logger, kegg_step_name, kegg_tmp_dir, config.error_file_path)
        write_done_file(logger, kegg_done_file_path)

    # 12_3.  add annotations to otrhogroups table
    step_number = '12_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_annotations'
    output_dir_path, tmp_dir = prepare_directories(logger, config.steps_results_dir, config.tmp_dir, step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Adding annotations to orthogroups...')
        start_time = time.time()

        final_orthologs_df = pd.read_csv(final_orthologs_table_file_path, dtype=str)

        if kegg_table_path.exists():
            kegg_table_df = pd.read_csv(kegg_table_path)[['OG_name', 'knum', 'knum_description']]
            final_orthologs_df = pd.merge(kegg_table_df, final_orthologs_df, on='OG_name')

        if cai_table_path.exists():
            cai_df = pd.read_csv(cai_table_path)[['OG_name', 'CAI_mean']]
            final_orthologs_df = pd.merge(cai_df, final_orthologs_df, on='OG_name')

        orthogroups_sizes_df = pd.read_csv(orthogroups_visualizations_dir_path / 'orthogroups_sizes.csv')
        final_orthologs_df = pd.merge(orthogroups_sizes_df, final_orthologs_df, on='OG_name')

        final_orthologs_table_annotated_path = output_dir_path / 'orthogroups_annotated.csv'
        final_orthologs_df.to_csv(final_orthologs_table_annotated_path, index=False)
        logger.info(f'Final orthologs table with annotations saved to {final_orthologs_table_annotated_path}')

        step_time = timedelta(seconds=int(time.time() - start_time))
        times_logger.info(f'Step {step_name} took {step_time}.')

        if not config.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, output_dir_path, config.final_output_dir)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    update_progressbar(config.progressbar_file_path, 'Analyze orthogroups codon bias')
    update_progressbar(config.progressbar_file_path, 'Annotate orthogroups with KO terms')


def run_main_pipeline(logger, times_logger, config):
    with open(config.genomes_names_path, 'r') as genomes_names_fp:
        genomes_names = genomes_names_fp.read().split('\n')

    if config.filter_out_plasmids or config.inputs_fasta_type == 'orfs':
        filtered_inputs_dir = step_1_fix_input_files(logger, times_logger, config)
        config.data_path = filtered_inputs_dir

    if config.step_to_complete == '1':
        logger.info("Step 1 completed.")
        return

    orfs_dir, translated_orfs_dir, all_orfs_fasta_path, all_proteins_fasta_path, orfs_coordinates_dir = \
        step_2_search_orfs(logger, times_logger, config)

    if config.step_to_complete == '2':
        logger.info("Step 2 completed.")
        return

    return

    if not config.only_calc_ogs:
        step_3_analyze_genome_completeness(logger, times_logger, config, translated_orfs_dir)

    if config.step_to_complete == '3':
        logger.info("Step 3 completed.")
        return

    final_orthogroups_file_path = step_5_orthogroups_inference(logger, times_logger, config, genomes_names,
                                                               translated_orfs_dir, all_proteins_fasta_path,
                                                               orfs_coordinates_dir)

    if config.step_to_complete == '5':
        logger.info("Step 5 completed.")
        return

    orthogroups_visualizations_dir_path = step_7_orthologs_table_variations(logger, times_logger, config,
                                                                            final_orthogroups_file_path, orfs_coordinates_dir)

    if config.step_to_complete == '7':
        logger.info("Step 7 completed.")
        return

    if config.only_calc_ogs:
        return

    ogs_dna_sequences_path, ogs_aa_sequences_path, ogs_aa_alignments_path, ogs_dna_alignments_path, \
        ogs_aa_consensus_dir_path = step_8_build_orthologous_groups_fastas(
        logger, times_logger, config, all_orfs_fasta_path, all_proteins_fasta_path, final_orthogroups_file_path)

    if config.step_to_complete == '8':
        logger.info("Step 8 completed.")
        return

    aligned_core_proteome_reduced_file_path, core_proteome_reduced_length = step_9_extract_core_genome_and_core_proteome(
        logger, times_logger, config, ogs_aa_alignments_path, ogs_dna_alignments_path)

    if config.step_to_complete == '9':
        logger.info("Step 9 completed.")
        return

    step_11_phylogeny(logger, times_logger, config, aligned_core_proteome_reduced_file_path, genomes_names,
                      core_proteome_reduced_length)

    if config.step_to_complete == '11':
        logger.info("Step 11 completed.")
        return

    step_12_orthogroups_annotations(logger, times_logger, config, orfs_dir, ogs_dna_sequences_path,
                                    ogs_aa_sequences_path, ogs_aa_consensus_dir_path, final_orthogroups_file_path,
                                    orthogroups_visualizations_dir_path)

    if config.step_to_complete == '12':
        logger.info("Step 12 completed.")
        return


def main():
    start_time = time.time()

    logger, times_logger, config = get_configuration()

    try:
        initialize_progressbar(config)
        prepare_and_verify_input_data(logger, config)
        update_progressbar(config.progressbar_file_path, 'Validate input files')

        run_main_pipeline(logger, times_logger, config)

        zip_results(logger, config)
        submit_clean_folders_job(logger, config)

        update_progressbar(config.progressbar_file_path, 'Finalize results')
        write_done_file(logger, config.done_files_dir / 'pipeline_finished_successfully.txt')
        state = State.Finished
    except subprocess.CalledProcessError as e:
        error_message = f'Error in command: "{e.cmd}": {e.stderr}'
        logger.exception(error_message)
        with open(config.error_file_path, 'a+') as f:
            f.write(error_message)
        state = State.Crashed
    except Exception as e:
        if not config.error_file_path.exists():
            with open(config.error_file_path, 'w') as f:
                f.write(traceback.format_exc())
        state = State.Crashed

    total_time = timedelta(seconds=int(time.time() - start_time))
    times_logger.info(f'Total pipeline time: {total_time}. Done')

    send_email_in_pipeline_end(logger, config, state)
    submit_clean_old_user_results_job(logger, config)


if __name__ == '__main__':
    print(f'sys.path is\n{sys.path}')
    main()
