import argparse
import logging
import math
import os
import shutil
import sys
from time import time
import traceback
import json
import subprocess

import matplotlib.pyplot as plt
import mmap
import numpy as np
import pandas as pd
import seaborn as sns

from auxiliaries.email_sender import send_email
from auxiliaries.file_writer import write_to_file
from auxiliaries.input_verifications import prepare_and_verify_input_data
from auxiliaries.pipeline_auxiliaries import measure_time, execute, wait_for_results, \
    prepare_directories, fail, submit_mini_batch, submit_batch, notify_admin, add_results_to_final_dir, remove_path,\
    str_to_bool
from auxiliaries.html_editor import edit_success_html, edit_failure_html, edit_progress
from auxiliaries import consts, cgi_consts
from flask import flask_interface_consts
from auxiliaries.logic_auxiliaries import mimic_prodigal_output, aggregate_ani_results, remove_bootstrap_values, \
    plot_genomes_histogram, update_progressbar, define_intervals
from flask.SharedConsts import USER_FILE_NAME_ZIP, USER_FILE_NAME_TAR

PIPELINE_STEPS = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--args_json_path', help='path to a json file that contains values for arguments which will '
                                                 'override the default values. Optional.')
    parser.add_argument('--run_dir', help='path to a directory where the pipeline will be run. Should contain a zip of'
                                          'the genomes. Mutually exclusive with --contigs_dir.')
    parser.add_argument('--contigs_dir',
                        help='path to a folder with the genomic sequences. This folder may be zipped, as well the files'
                             ' in it. Mutually exclusive with --run_dir.')
    parser.add_argument('--output_dir', help='relative path of directory where the output files will be written to',
                        default='outputs')
    parser.add_argument('--email', help='A notification will be sent once the pipeline is done',
                        default=flask_interface_consts.OWNER_EMAIL)
    parser.add_argument('--identity_cutoff', default=40,
                        help='minimum required percent of identity level (lower values will be filtered out)')
    parser.add_argument('--coverage_cutoff', default=70,
                        help='minimum required coverage percent of homology region between genes (lower values will be filtered out)')
    parser.add_argument('--e_value_cutoff', default=0.01,
                        help='maxmimum permitted e-value (0 <= e_value_cutoff <= 1; higher values will be filtered'
                             ' out).')
    parser.add_argument('--core_minimal_percentage', default=100.0,
                        help='the minimum required percent of gene members that is needed to be considered a core gene.'
                             ' For example: (1) 100 means that for a gene to be considered core, all strains should '
                             'have a member in the group.\n(2) 50 means that for a gene to be considered core, at least'
                             ' half of the strains should have a member in the group.\n(3) 0 means that every gene '
                             'should be considered as a core gene.')
    parser.add_argument('--bootstrap', default='no', choices=['yes', 'no'],
                        help='whether or not to apply bootstrap procedure over the reconstructed species tree.')
    parser.add_argument('--outgroup', default=None,
                        help='The species name used to to root the phylogenetic tree, or None to leave unrooted.')
    parser.add_argument('--filter_out_plasmids', action='store_true',
                        help='whether or not to filter out plasmids from the input files')
    parser.add_argument('--inputs_fasta_type', choices=['genomes', 'orfs'], default='genomes',
                        help='whether the input files are fastas of orfs and therefore Prodigal is skipped, '
                             'or genomes assemblies and then the first step is ORFs extraction using Prodigal')
    parser.add_argument('--add_orphan_genes_to_ogs', action='store_true',
                        help='whether orphan genes should be considered as OGs')
    parser.add_argument('--qfo_benchmark', action='store_true',
                        help='whether the input files are annotated genomes in the QfO benchmark format')
    # choices=['pupkoweb', 'pupkowebr', 'pupkolab', 'pupkolabr', 'pupkotmp', 'pupkotmpr', 'itaym', 'lilach',
    # 'bioseq', 'bental', 'oren.q', 'bioseq20.q'])
    parser.add_argument('-q', '--queue_name', help='The queue to which the job(s) will be submitted to',
                        default=consts.DEFAULT_SLURM_PARTITION)
    parser.add_argument('--account_name', help='The slurm account to submit jobs to',
                        default=consts.DEFAULT_SLURM_ACCOUNT)
    parser.add_argument('--step_to_complete', help='The final step to execute', default=None,
                        choices=[*PIPELINE_STEPS, None])
    parser.add_argument('--only_calc_ogs', help='Do only the necessary steps to calculate OGs', action='store_true')
    parser.add_argument('--zip_results_in_partial_pipeline', help='Zip results also when the pipeline is partially run', action='store_true')
    parser.add_argument('--bypass_number_of_genomes_limit', help='Bypass the limit on number of genomes',
                        action='store_true')
    parser.add_argument('--optimize_orthogroups_inference', help='Optimize the orthogroups inference using heuristics',
                        action='store_true')
    parser.add_argument('--num_of_clusters_in_orthogroup_inference', help='Number of clusters in the optimization of '
                                                                          'the orthogroups inference. Relevant only if '
                                                                          'optimize_orthogroups_inference is True',
                        default=5)
    parser.add_argument('--run_optimized_mmseqs', help='Optimize the mmseqs run',
                        action='store_true')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

    # parser.add_argument('--promoters_length', default=300,
    #                     help='How many basepairs upstream to the ATG will be used for the sweeps analysis',
    #                     type=lambda x: int(x) if int(x) >= 0 else parser.error(
    #                         f'Minimal number of upstream basepairs should be non-negative!'))
    # parser.add_argument('--minimal_number_of_sequences_allowed_for_sweeps_analysis', default=21,
    #                     help='MSAs with fewer sequences will be ignore in the sweeps analysis.',
    #                     type=lambda x: int(x) if int(x) > 4 else parser.error(
    #                         f'Minimal number of sequences required for sweeps analysis is 5!'))

    args = parser.parse_args()

    # Override arguments with args_json_path content
    if args.args_json_path:
        with open(args.args_json_path, 'r') as args_json_file:
            args_json = json.load(args_json_file)
            args.__dict__.update(args_json)

    return args


def prepare_pipeline_framework(args):
    if (args.run_dir and args.contigs_dir) or (not args.run_dir and not args.contigs_dir):
        raise ValueError('Either run_dir or contigs_dir should be provided, but not both.')

    if args.run_dir:
        meta_output_dir = args.run_dir
        if os.path.exists(os.path.join(meta_output_dir, USER_FILE_NAME_ZIP)):
            args.contigs_dir = os.path.join(meta_output_dir, USER_FILE_NAME_ZIP)
        elif os.path.exists(os.path.join(meta_output_dir, USER_FILE_NAME_TAR)):
            args.contigs_dir = os.path.join(meta_output_dir, USER_FILE_NAME_TAR)
        else:
            raise ValueError(f'No genomes zip or tar file found in {meta_output_dir}')
    else:  # args.contigs_dir was provided
        if os.path.exists(args.contigs_dir):
            args.contigs_dir = args.contigs_dir.rstrip('/')
            meta_output_dir = os.path.dirname(args.contigs_dir)
        else:
            raise ValueError(f'contigs_dir argument {args.contigs_dir} does not exist!')

    output_dir = os.path.join(meta_output_dir, args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(format=consts.LOG_MESSAGE_FORMAT,
                        level=level)
    logger = logging.getLogger('main')
    times_logger = logging.getLogger('times')

    if consts.LOG_IN_SEPARATE_FILES:
        formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)

        main_file_handler = logging.FileHandler(os.path.join(output_dir, 'main_log.txt'))
        main_file_handler.setFormatter(formatter)
        logger.addHandler(main_file_handler)

        times_file_handler = logging.FileHandler(os.path.join(output_dir, 'times_log.txt'))
        times_file_handler.setFormatter(formatter)
        times_logger.addHandler(times_file_handler)

    logger.info(args)
    logger.info(f'meta_output_dir is: {meta_output_dir}')
    logger.info(f'Created output_dir in: {output_dir}')

    error_file_path = os.path.join(meta_output_dir, flask_interface_consts.ERROR_FILE_PATH)
    progressbar_file_path = os.path.join(meta_output_dir, flask_interface_consts.PROGRESSBAR_FILE_PATH)

    run_number = os.path.join(os.path.split(meta_output_dir)[1])
    logger.info(f'run_number is {run_number}')

    if not consts.IGNORE_HTML:
        output_html_path = os.path.join(meta_output_dir, cgi_consts.RESULT_WEBPAGE_NAME)
        logger.info(f'output_html_path is {output_html_path}')

        with open(output_html_path) as f:
            html_content = f.read()
        html_content = html_content.replace('QUEUED', 'RUNNING')
        if 'progress-bar-striped active' not in html_content:
            html_content = html_content.replace('progress-bar-striped', 'progress-bar-striped active')
        with open(output_html_path, 'w') as f:
            f.write(html_content)

        output_url = os.path.join(cgi_consts.WEBSERVER_RESULTS_URL, run_number, cgi_consts.RESULT_WEBPAGE_NAME)
        logger.info(f'output_url is {output_url}')

        meta_output_url = os.path.join(cgi_consts.WEBSERVER_RESULTS_URL, run_number)
    else:
        output_html_path = ''
        output_url = ''
        meta_output_url = ''

    tmp_dir = os.path.join(output_dir, 'tmp')
    logger.info(f'Creating tmp_dir in: {tmp_dir}')
    os.makedirs(tmp_dir, exist_ok=True)

    done_files_dir = os.path.join(output_dir, 'done')
    logger.info(f'Creating done_files_dir in: {done_files_dir}')
    os.makedirs(done_files_dir, exist_ok=True)

    steps_results_dir = os.path.join(output_dir, 'steps_results')
    logger.info(f'Creating results_dir in: {steps_results_dir}')
    os.makedirs(steps_results_dir, exist_ok=True)

    data_path = os.path.join(output_dir, 'inputs')
    logger.info(f'Creating data_path is: {data_path}')
    os.makedirs(data_path, exist_ok=True)

    return logger, times_logger, meta_output_dir, error_file_path, progressbar_file_path, run_number, output_html_path, \
        output_url, meta_output_url, output_dir, tmp_dir, done_files_dir, steps_results_dir, data_path


def validate_arguments(args):
    if type(args.bootstrap) == str:
        args.bootstrap = str_to_bool(args.bootstrap)
    if type(args.filter_out_plasmids) == str:
        args.filter_out_plasmids = str_to_bool(args.filter_out_plasmids)
    if type(args.add_orphan_genes_to_ogs) == str:
        args.add_orphan_genes_to_ogs = str_to_bool(args.add_orphan_genes_to_ogs)
    if type(args.optimize_orthogroups_inference) == str:
        args.optimize_orthogroups_inference = str_to_bool(args.optimize_orthogroups_inference)

    if args.outgroup == "No outgroup":
        args.outgroup = None

    args.identity_cutoff = float(args.identity_cutoff)
    if args.identity_cutoff < 0 or args.identity_cutoff > 100:
        raise ValueError(f'identity_cutoff argument {args.identity_cutoff} has invalid value')

    args.coverage_cutoff = float(args.coverage_cutoff)
    if args.coverage_cutoff < 0 or args.coverage_cutoff > 100:
        raise ValueError(f'coverage_cutoff argument {args.coverage_cutoff} has invalid value')

    args.e_value_cutoff = float(args.e_value_cutoff)
    if args.e_value_cutoff < 0 or args.e_value_cutoff > 1:
        raise ValueError(f'e_value_cutoff argument {args.e_value_cutoff} has invalid value')

    args.core_minimal_percentage = float(args.core_minimal_percentage)
    if args.core_minimal_percentage < 0 or args.core_minimal_percentage > 100:
        raise ValueError(f'core_minimal_percentage argument {args.core_minimal_percentage} has invalid value')

    args.num_of_clusters_in_orthogroup_inference = int(args.num_of_clusters_in_orthogroup_inference)
    if args.num_of_clusters_in_orthogroup_inference <= 1:
        raise ValueError(f'num_of_clusters_in_orthogroup_inference argument {args.num_of_clusters_in_orthogroup_inference} has invalid value')

    if not args.optimize_orthogroups_inference:
        args.num_of_clusters_in_orthogroup_inference = 1


def initialize_progressbar(args, progressbar_file_path):
    if args.only_calc_ogs:
        steps = consts.ONLY_CALC_OGS_TABLE_STEPS_NAMES_FOR_PROGRESS_BAR
    else:
        steps = consts.FULL_STEPS_NAMES_FOR_PROGRESS_BAR

    if not args.filter_out_plasmids:
        steps.remove('Filter out plasmids')

    if args.inputs_fasta_type == 'orfs':
        for i in range(len(steps)):
            if steps[i] == 'Predict and translate ORFs':
                steps[i] = 'Translate ORFs'
                break

    if args.core_minimal_percentage != 100:
        steps.remove('Calculate genomes numeric representation')

    df = pd.DataFrame({'Step': steps, 'Finished': [False] * len(steps)})
    df.to_csv(progressbar_file_path, index=False)


def step_0_fix_input_files(args, logger, times_logger, error_file_path, output_dir, tmp_dir, done_files_dir,
                           data_path):
    # 0.	drop_plasmids.py
    step_number = '00'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_fix_input_files'
    script_path = os.path.join(consts.SRC_DIR, 'steps/drop_plasmids_and_fix_frames.py')
    filtered_inputs_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Filtering plasmids out...')
        all_cmds_params = []
        for fasta_file in os.listdir(data_path):
            single_cmd_params = [os.path.join(data_path, fasta_file),
                                 os.path.join(filtered_inputs_dir, fasta_file)]

            if args.filter_out_plasmids:
                single_cmd_params.append('--drop_plasmids')
            if args.inputs_fasta_type == 'orfs':
                single_cmd_params.append('--fix_frames')

            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='drop_plasmids',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return filtered_inputs_dir


def step_1_calculate_ani(args, logger, times_logger, error_file_path,  output_dir, tmp_dir, final_output_dir,
                         done_files_dir, data_path):
    # 1.   ani.py
    step_number = '01'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_ani'
    script_path = os.path.join(consts.SRC_DIR, 'steps/ani.py')
    ani_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Calculating ANI values...')
        ani_tmp_files = os.path.join(pipeline_step_tmp_dir, 'temp_results')
        os.makedirs(ani_tmp_files, exist_ok=True)

        genomes_paths = [os.path.join(data_path, genome_file_name) for genome_file_name in os.listdir(data_path)]
        genomes_list_path = os.path.join(ani_tmp_files, 'genomes_list.txt')
        with open(genomes_list_path, 'w') as genomes_list_file:
            genomes_list_file.write('\n'.join(genomes_paths))

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for fasta_file in os.listdir(data_path):
            single_cmd_params = [os.path.join(data_path, fasta_file), genomes_list_path, ani_tmp_files]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='calculate_ani',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   memory=consts.ANI_REQUIRED_MEMORY_GB)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        # Aggregate ANI results
        aggregate_ani_results(ani_tmp_files, ani_output_dir)

        add_results_to_final_dir(logger, ani_output_dir, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def step_2_search_orfs(args, logger, times_logger, error_file_path,  output_dir, tmp_dir, final_output_dir,
                       done_files_dir, data_path):
    # 2a.	search_orfs.py
    # Input: (1) an input path for a fasta file with contigs/full genome (2) an output file path
    #                (with a suffix as follows: i_genes.fasta. especially relevant for the wrapper).
    # Output: a fasta file where under each header, there’s a single gene.
    # Can be parallelized on cluster
    # Prodigal ignores newlines in the middle of a sequence, namely, >bac1\nAAA\nAA\nTTT >bac1\nAAAAATTT will be
    # analyzed identically.
    step_number = '02a'  # using orfs_step_number for downstream analysis (extract promoter and gene sequences)
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orfs'
    script_path = os.path.join(consts.SRC_DIR, 'steps/search_orfs.py')
    orfs_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        if args.inputs_fasta_type == 'genomes':
            logger.info('Extracting ORFs...')
            all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
            for fasta_file in os.listdir(data_path):
                single_cmd_params = [f'"{os.path.join(data_path, fasta_file)}"',
                                     orfs_dir,
                                     step_name]
                all_cmds_params.append(single_cmd_params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                       num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                       job_name_suffix='search_orfs',
                                                       queue_name=args.queue_name,
                                                       account_name=args.account_name)

            wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                             num_of_batches, error_file_path)
        else:  # inputs are orfs
            logger.info(f'Inputs are already annotated genomes. Skipping step {step_name}.')
            shutil.copytree(data_path, orfs_dir, dirs_exist_ok=True)
            mimic_prodigal_output(orfs_dir, step_name)

        add_results_to_final_dir(logger, orfs_dir, final_output_dir, keep_in_source_dir=True)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # make sure that all ORF files contain something....
    missing_orfs = 0
    error_msg = ''
    for file in sorted(os.listdir(orfs_dir)):
        if '02_ORFs' not in file:  # Ignore system files that are automatically created sometimes
            continue
        try:
            with open(os.path.join(orfs_dir, file), 'rb', 0) as orf_f, mmap.mmap(orf_f.fileno(), 0,
                                                                                 access=mmap.ACCESS_READ) as s:
                if s.find(b'>') > -1:
                    continue
        except:
            if not missing_orfs:
                # add msg prefix when a genome without orfs detected for the first time
                error_msg = f'{flask_interface_consts.WEBSERVER_NAME} could not detect any ORFs in:\n<br> '
                missing_orfs = 1
            # add genome-without-orfs name
            error_msg += f'{os.path.splitext(file)[0]}\n<br>'

    if missing_orfs:
        error_msg += f'\n<br>Please remove the abovementioned from the dataset and re-submit your job. In general, ' \
                     f'it is recommended to use genomic files that contain at least 20K base pairs (each).'
        fail(logger, error_msg, error_file_path)

    # 2b.	extract_orfs_statistics.py
    # Input: (1) A path to ORFs file (2) An output dir of orfs statistics
    # Output: write the number of ORFs and GC content to the output files (respectively)
    # Can be parallelized on cluster
    step_number = '02b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orfs_statistics'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_orfs_statistics.py')
    orfs_statistics_dir, orfs_statistics_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Collecting orfs counts...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for file in os.listdir(orfs_dir):
            orf_path = os.path.join(orfs_dir, file)
            single_cmd_params = [orf_path, orfs_statistics_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_expected_orfs_results, example_cmd = submit_batch(logger, script_path, all_cmds_params,
                                                                 orfs_statistics_tmp_dir, error_file_path,
                                                                 num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                                 job_name_suffix='orfs_statistics',
                                                                 queue_name=args.queue_name,
                                                                 account_name=args.account_name,)

        wait_for_results(logger, times_logger, step_name, orfs_statistics_tmp_dir,
                         num_of_expected_orfs_results, error_file_path=error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 2c.	plot_orfs_statistics
    step_number = '02c'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orfs_plots'
    orfs_plots_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Concatenating orfs counts and gc contents...')

        gc_content = {}
        orf_count = {}
        for file_name in os.listdir(orfs_statistics_dir):
            strain_name = os.path.splitext(file_name)[0]
            with open(os.path.join(orfs_statistics_dir, file_name), 'r') as fp:
                orfs_statistics = json.load(fp)

            gc_content[strain_name] = orfs_statistics['gc_content']
            orf_count[strain_name] = orfs_statistics['orfs_count']

        plot_genomes_histogram(orf_count, orfs_plots_path, 'orfs_counts', 'ORFs count', 'ORFs Count per genome')
        plot_genomes_histogram(gc_content, orfs_plots_path, 'orfs_gc_content', 'GC Content (of ORFs)', 'GC Content per genome')

        add_results_to_final_dir(logger, orfs_plots_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 2d.  translate_fna_to_faa.py - of ORFs files
    # Input: path to fna file and a faa file
    # Output: translate the fna to protein and write to the faa file
    step_number = '02d'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_translated_orfs'
    script_path = os.path.join(consts.SRC_DIR, 'steps/translate_fna_to_faa.py')
    translated_orfs_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Translating ORFs sequences...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for fasta_file in os.listdir(orfs_dir):
            file_path = os.path.join(orfs_dir, fasta_file)
            translated_file_name = f'{os.path.splitext(fasta_file)[0]}.{step_name}'
            output_path = os.path.join(translated_orfs_dir_path, translated_file_name)

            single_cmd_params = [file_path, output_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='orfs_dna_translation',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        add_results_to_final_dir(logger, translated_orfs_dir_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 2e.  Aggregate all protein fasta files to one file
    step_number = '02e'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_concat_translated_orfs'
    concat_translated_orfs_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    all_proteins_fasta_path = os.path.join(concat_translated_orfs_dir_path, 'all_proteomes.faa')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Concatenating translated ORFs sequences...')

        cmd = f"cat {os.path.join(translated_orfs_dir_path, '*')} > {all_proteins_fasta_path}"
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orfs_dir, translated_orfs_dir_path, all_proteins_fasta_path


def step_3_analyze_genome_completeness(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                       final_output_dir, done_files_dir, translated_orfs_dir_path):
    # 3.  assessing_genomes_completeness.py
    # Input: path to faa file
    # Output: calculate proteome completeness based on a dataset of profile-HMMs that represent core bacterial genes.
    step_number = '03'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_genomes_completeness'
    script_path = os.path.join(consts.SRC_DIR, 'steps/assessing_genome_completeness.py')
    genome_completeness_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    genomes_output_dir_path = os.path.join(genome_completeness_dir_path, 'individual_proteomes_outputs')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Calculating genomes completeness...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        os.makedirs(genomes_output_dir_path, exist_ok=True)
        for fasta_file in os.listdir(translated_orfs_dir_path):
            file_path = os.path.join(translated_orfs_dir_path, fasta_file)
            single_cmd_params = [file_path, genomes_output_dir_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='genomes_completeness',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        # Aggregate results of step
        genomes_completeness_scores = {}
        for genome_name in os.listdir(genomes_output_dir_path):
            genome_score_path = os.path.join(genomes_output_dir_path, genome_name, 'result.txt')
            with open(genome_score_path, 'r') as fp:
                genomes_completeness_scores[genome_name] = float(fp.read().strip())

        plot_genomes_histogram(genomes_completeness_scores, genome_completeness_dir_path, 'genomes_completeness',
                               'Genome BUSCO completeness score', 'Genome BUSCO completeness score per genome')

        # comment the next line if you don't wish to delete hmmer results
        shutil.rmtree(genomes_output_dir_path)

        add_results_to_final_dir(logger, genome_completeness_dir_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def step_4_cluster_proteomes(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                            done_files_dir, all_proteins_fasta_path):
    # 4.	cluster_proteomes.py
    step_number = '04'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_cluster_proteomes'
    script_path = os.path.join(consts.SRC_DIR, 'steps/cluster_proteomes.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    clusters_file_path = os.path.join(pipeline_step_output_dir, 'clusters.tsv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        params = [all_proteins_fasta_path,
                  pipeline_step_output_dir,
                  clusters_file_path,
                  consts.MMSEQS_CLUSTER_MIN_SEQ_ID,
                  consts.MMSEQS_CLUSTER_MIN_COVERAGE,
                  consts.MMSEQS_CLUSTER_NUM_OF_CORES,
                  f'--num_of_clusters_in_orthogroup_inference {args.num_of_clusters_in_orthogroup_inference}'
                  ]

        if args.optimize_orthogroups_inference:
            params.append('--do_cluster')

        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path, args.queue_name, args.account_name,
                          job_name='cluster_proteomes', num_of_cpus=consts.MMSEQS_CLUSTER_NUM_OF_CORES, memory=consts.MMSEQS_REQUIRED_MEMORY_GB)
        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return pipeline_step_output_dir


def step_5_infer_orthogroups(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                             done_files_dir, clusters_dir, translated_orfs_dir, all_proteins_path, genomes_names_path):
    # 50.  filter_fasta_files_by_cluster.py
    step_number = '050'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_filter_fastas_by_cluster'
    script_path = os.path.join(consts.SRC_DIR, 'steps/filter_fasta_file.py')
    clusters_fasta_files, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs

        for filename in os.listdir(clusters_dir):
            if 'cluster_' not in filename:
                continue
            cluster_genes_path = os.path.join(clusters_dir, filename)
            cluster_id = os.path.splitext(filename)[0].split('_')[1]
            cluster_fasta_files = os.path.join(clusters_fasta_files, f'cluster_{cluster_id}_proteomes')
            os.makedirs(cluster_fasta_files, exist_ok=True)

            if not args.run_optimized_mmseqs:
                for protein_fasta_file in os.listdir(translated_orfs_dir):
                    params = [
                        cluster_genes_path,
                        os.path.join(translated_orfs_dir, protein_fasta_file),
                        os.path.join(cluster_fasta_files, protein_fasta_file)
                    ]
                    all_cmds_params.append(params)
            else:
                params = [cluster_genes_path, all_proteins_path, os.path.join(cluster_fasta_files, 'all_proteomes.faa')]
                all_cmds_params.append(params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='filter_fastas_by_clusters',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


    # 5.  infer_orthogroups.py
    step_number = '05'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_infer_orthogroups'
    script_path = os.path.join(consts.SRC_DIR, 'steps/infer_orthogroups.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    pipeline_step_done_dir = os.path.join(done_files_dir, step_name)
    os.makedirs(pipeline_step_done_dir, exist_ok=True)
    orthogroups_file_path = os.path.join(pipeline_step_output_dir, 'orthogroups.csv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs

        for cluster_fastas_dir_name in os.listdir(clusters_fasta_files):
            cluster_fastas_dir_path = os.path.join(clusters_fasta_files, cluster_fastas_dir_name)
            cluster_id = cluster_fastas_dir_name.split('_')[1]

            cluster_output_dir = os.path.join(pipeline_step_output_dir, str(cluster_id))
            os.makedirs(cluster_output_dir, exist_ok=True)
            cluster_tmp_dir = os.path.join(pipeline_step_tmp_dir, str(cluster_id))
            os.makedirs(cluster_tmp_dir, exist_ok=True)
            cluster_done_dir = os.path.join(pipeline_step_done_dir, str(cluster_id))
            os.makedirs(cluster_done_dir, exist_ok=True)

            params = [step_number, cluster_output_dir, cluster_tmp_dir, cluster_done_dir, cluster_fastas_dir_path,
                      os.path.join(cluster_fastas_dir_path, 'all_proteomes.faa'), genomes_names_path, args.queue_name,
                      args.account_name, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff,
                      args.num_of_clusters_in_orthogroup_inference]

            if args.run_optimized_mmseqs:
                params.append('--run_optimized_mmseqs')

            all_cmds_params.append(params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='infer_orthogroups',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   memory=consts.INFER_ORTHOGROUPS_MEMORY_GB,
                                                   time_in_hours=96)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, recursive_step=True)

        # Aggregate OG tables of all clusters
        all_og_tables = []
        for cluster_dir_name in os.listdir(pipeline_step_output_dir):
            cluster_og_table = os.path.join(pipeline_step_output_dir, cluster_dir_name, '05j_verified_table', 'verified_orthologs_table.csv')
            cluster_og_df = pd.read_csv(cluster_og_table)
            all_og_tables.append(cluster_og_df)
        unified_og_table = pd.concat(all_og_tables, axis=0, ignore_index=True)
        unified_og_table['OG_name'] = [f'OG_{i}' for i in range(len(unified_og_table.index))]
        unified_og_table.to_csv(orthogroups_file_path, index=False)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthogroups_file_path


def step_6_extract_orphan_genes(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                                done_files_dir, orfs_dir, orthologs_table_file_path):
    # 6.   extract_orphan_genes.py
    step_number = '06'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orphan_genes'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_orphan_genes.py')
    orphan_genes_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orphan_genes_internal_dir = os.path.join(orphan_genes_dir, 'orphans_lists_per_genome')
    os.makedirs(orphan_genes_internal_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting orphan genes...')
        all_cmds_params = []

        for orf_file in os.listdir(orfs_dir):
            single_cmd_params = [os.path.join(orfs_dir, orf_file), orthologs_table_file_path, orphan_genes_internal_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='extract_orphans',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        all_stat_dfs = []
        for file_name in os.listdir(orphan_genes_internal_dir):
            if 'orphans_stats.csv' not in file_name:
                continue
            df = pd.read_csv(os.path.join(orphan_genes_internal_dir, file_name), index_col=0)
            all_stat_dfs.append(df)

        combined_df = pd.concat(all_stat_dfs)
        combined_df.to_csv(os.path.join(orphan_genes_dir, 'orphans_genes_stats.csv'))

        number_of_orphans_per_file = combined_df['Total orphans count'].to_dict()
        plot_genomes_histogram(number_of_orphans_per_file, orphan_genes_dir, 'orphan_genes_count', 'Orphan genes count',
                               'Orphan genes count per Genome')

        add_results_to_final_dir(logger, orphan_genes_dir, final_output_dir, keep_in_source_dir=True)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orphan_genes_internal_dir


def step_7_orthologs_table_variations(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                      final_output_dir, done_files_dir, orthologs_table_file_path, orphan_genes_dir):
    # 7a.	orthologs_table_variations.py
    # Input: (1) path to the orthologs table (2) path to the orphan genes directory.
    # Output: build variations of the orthologs table.
    step_number = '07a'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups'
    script_path = os.path.join(consts.SRC_DIR, 'steps/orthologs_table_variations.py')
    final_orthologs_table_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    final_orthologs_table_file_path = os.path.join(final_orthologs_table_dir_path, 'orthogroups.csv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing final orthologs table...')
        params = [orthologs_table_file_path,
                  final_orthologs_table_file_path
                  ]
        if args.qfo_benchmark:
            params += ['--qfo_benchmark']
        if args.add_orphan_genes_to_ogs:
            params += [f'--orphan_genes_dir {orphan_genes_dir}']
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path, args.queue_name, args.account_name,
                          job_name='final_ortholog_groups')
        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        add_results_to_final_dir(logger, final_orthologs_table_dir_path, final_output_dir,
                                 keep_in_source_dir=True)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 7b.	extract_groups_sizes_frequency
    step_number = '07b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_sizes'
    group_sizes_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Collecting sizes...')

        final_orthologs_table_df = pd.read_csv(final_orthologs_table_file_path, index_col='OG_name')
        group_sizes = final_orthologs_table_df.apply(lambda row: row.count(), axis=1)
        group_sizes.name = 'OG size (number of genomes)'
        group_sizes.to_csv(os.path.join(group_sizes_path, 'groups_sizes.csv'))

        sns.histplot(x=group_sizes, discrete=True)
        if len(np.unique(group_sizes)) < 10:
            plt.xticks(np.unique(group_sizes))
        plt.title('Orthologous groups sizes distribution', fontsize=20, loc='center', wrap=True)
        plt.xlabel('OG size (number of genomes)', fontsize=15)
        plt.ylabel('Count', fontsize=15)
        plt.tight_layout()
        plt.savefig(os.path.join(group_sizes_path, 'groups_sizes.png'), dpi=600)
        plt.clf()

        add_results_to_final_dir(logger, group_sizes_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return final_orthologs_table_file_path


def step_8_build_orthologous_groups_fastas(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                           final_output_dir, done_files_dir, orfs_dir, final_orthologs_table_file_path):
    # 8a.	extract_orfs.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    step_number = '08a'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_dna'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_orfs.py')
    orthologs_dna_sequences_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir,
                                                                                  step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        with open(final_orthologs_table_file_path, 'r') as fp:
            number_of_ogs = sum(1 for _ in fp) - 1

        ogs_to_process_per_job = math.ceil(number_of_ogs / consts.MAX_PARALLEL_JOBS)
        lines_intervals = define_intervals(0, number_of_ogs, ogs_to_process_per_job)
        for (start_index, end_index_exclusive) in lines_intervals:
            single_cmd_params = [orfs_dir,
                                 final_orthologs_table_file_path,
                                 start_index,
                                 end_index_exclusive,
                                 orthologs_dna_sequences_dir_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='orfs_extraction',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        add_results_to_final_dir(logger, orthologs_dna_sequences_dir_path, final_output_dir,
                                 keep_in_source_dir=True)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 8b.  translate_fna_to_faa.py
    # Input: path to fna file and an faa file
    # Output: translate the fna to protein and write to the faa file
    # Can be parallelized on cluster
    step_number = '08b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_aa'
    script_path = os.path.join(consts.SRC_DIR, 'steps/translate_fna_to_faa.py')
    orthologs_aa_sequences_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Translating orthologs groups sequences...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for fasta_file in os.listdir(orthologs_dna_sequences_dir_path):
            file_path = os.path.join(orthologs_dna_sequences_dir_path, fasta_file)
            output_path = os.path.join(orthologs_aa_sequences_dir_path, fasta_file.replace('_dna.fas', '_aa.fas'))

            single_cmd_params = [file_path, output_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='dna_translation',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        add_results_to_final_dir(logger, orthologs_aa_sequences_dir_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 8c.	align_orthologs_group.py
    # Input: (1) A path to an unaligned amino acid sequences file (2) An output file path
    # Output: aligned sequences
    # Can be parallelized on cluster
    step_number = '08c'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_aa_msa'
    script_path = os.path.join(consts.SRC_DIR, 'steps/align_orthologs_group.py')
    aa_alignments_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Aligning orthologs groups...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for og_file in os.listdir(orthologs_aa_sequences_dir_path):
            og_path = os.path.join(orthologs_aa_sequences_dir_path, og_file)
            og_file_prefix = os.path.splitext(og_file)[0]
            alignment_path = os.path.join(aa_alignments_path, f'{og_file_prefix}_mafft.fas')

            single_cmd_params = [og_path, alignment_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params,
                                                   pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='genes_alignment',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        add_results_to_final_dir(logger, aa_alignments_path, final_output_dir,
                                 keep_in_source_dir=True)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 8d.	induce_dna_msa_by_aa_msa.py
    # Input: (1) An aa alignment (2) An unaligned dna file (3) An output file path
    # Output: write the codon alignment induced by the aa alignment to the output file
    # Can be parallelized on cluster
    step_number = '08d'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_induced_dna_msa_by_aa_msa'
    script_path = os.path.join(consts.SRC_DIR, 'steps/induce_dna_msa_by_aa_msa.py')
    dna_alignments_path, induced_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info(f'Inducing dna alignments...\n(from {aa_alignments_path})')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for og_file in os.listdir(aa_alignments_path):
            aa_alignment_path = os.path.join(aa_alignments_path, og_file)
            dna_unaligned_path = os.path.join(orthologs_dna_sequences_dir_path, og_file.replace('aa_mafft', 'dna'))
            dna_induced_alignment_path = os.path.join(dna_alignments_path, og_file.replace('aa_mafft', 'dna_induced'))

            single_cmd_params = [aa_alignment_path, dna_unaligned_path, dna_induced_alignment_path]
            all_cmds_params.append(single_cmd_params)

        num_of_expected_induced_results, example_cmd = submit_batch(logger, script_path, all_cmds_params,
                                                                    induced_tmp_dir, error_file_path,
                                                                    num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                                    job_name_suffix='induce_msa',
                                                                    queue_name=args.queue_name,
                                                                    account_name=args.account_name)

        wait_for_results(logger, times_logger, step_name, induced_tmp_dir,
                         num_of_expected_results=num_of_expected_induced_results, error_file_path=error_file_path)

        add_results_to_final_dir(logger, dna_alignments_path, final_output_dir, keep_in_source_dir=True)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthologs_dna_sequences_dir_path, orthologs_aa_sequences_dir_path, aa_alignments_path, dna_alignments_path


def step_9_extract_core_genome_and_core_proteome(args, logger, times_logger, error_file_path,
                                                  output_dir, tmp_dir, final_output_dir, done_files_dir, num_of_strains,
                                                  strains_names_path, aa_alignments_path, dna_alignments_path):
    # 9a.	extract aligned_core_proteome.py
    step_number = '09a'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_aligned_core_proteome'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_core_genome.py')
    aligned_core_proteome_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    aligned_core_proteome_file_path = os.path.join(aligned_core_proteome_path, 'aligned_core_proteome.fas')
    core_proteome_ogs_names_file_path = os.path.join(aligned_core_proteome_path, 'core_ortholog_groups_names.txt')
    core_proteome_length_file_path = os.path.join(aligned_core_proteome_path, 'core_length.txt')
    number_of_core_proteome_members_file_path = os.path.join(aligned_core_proteome_path, 'number_of_core_genes.txt')

    if not os.path.exists(done_file_path):
        logger.info('Extracting aligned core proteome...')

        params = [aa_alignments_path, num_of_strains, strains_names_path,
                  aligned_core_proteome_file_path,
                  core_proteome_ogs_names_file_path,
                  core_proteome_length_file_path,
                  number_of_core_proteome_members_file_path,
                  f'--core_minimal_percentage {args.core_minimal_percentage}']  # how many members induce a core group?
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path,
                          args.queue_name, args.account_name, job_name='core_proteome')

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        add_results_to_final_dir(logger, aligned_core_proteome_path, final_output_dir,
                                 keep_in_source_dir=True)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    with open(core_proteome_length_file_path, 'r') as fp:
        core_proteome_length = int(fp.read().strip())


    # 9b.      extract_aligned_core_genome
    step_number = '09b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_aligned_core_genome'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_core_genome.py')
    aligned_core_genome_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    aligned_core_genome_file_path = os.path.join(aligned_core_genome_path, 'aligned_core_genome.fas')
    core_genome_ogs_names_file_path = os.path.join(aligned_core_genome_path, 'core_ortholog_groups_names.txt')
    core_genome_length_file_path = os.path.join(aligned_core_genome_path, 'core_length.txt')
    number_of_core_genome_members_file_path = os.path.join(aligned_core_genome_path, 'number_of_core_genes.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting aligned core genome...')

        params = [dna_alignments_path, num_of_strains, strains_names_path,
                  aligned_core_genome_file_path,
                  core_genome_ogs_names_file_path,
                  core_genome_length_file_path,
                  number_of_core_genome_members_file_path,
                  f'--core_minimal_percentage {args.core_minimal_percentage}']  # how many members induce a core group?
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path,
                          args.queue_name, args.account_name, job_name='core_genome')

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        add_results_to_final_dir(logger, aligned_core_genome_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return aligned_core_proteome_file_path, core_proteome_length


def step_10_genome_numeric_representation(args, logger, times_logger, error_file_path, output_dir,
                                         tmp_dir, final_output_dir, done_files_dir, orfs_dir,
                                         final_orthologs_table_file_path):
    # 10.	genome_numeric_representation.py
    step_number = '10'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_genome_numeric_representation'
    script_path = os.path.join(consts.SRC_DIR, 'steps/genome_numeric_representation.py')
    numeric_representation_output_dir, numeric_representation_tmp_dir = prepare_directories(logger, output_dir,
                                                                                            tmp_dir, step_name)
    core_genome_numeric_representation_file_path = os.path.join(numeric_representation_output_dir,
                                                                'core_genome_numeric_representation.txt')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        params = [final_orthologs_table_file_path,
                  orfs_dir,
                  core_genome_numeric_representation_file_path,
                  numeric_representation_tmp_dir
                  ]
        submit_mini_batch(logger, script_path, [params], numeric_representation_tmp_dir, error_file_path,
                          args.queue_name, args.account_name, job_name='numeric_representation')

        wait_for_results(logger, times_logger, step_name, numeric_representation_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        add_results_to_final_dir(logger, numeric_representation_output_dir, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def step_11_phylogeny(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                      done_files_dir, aligned_core_proteome_file_path, genomes_names_path, number_of_genomes, core_proteome_length):
    # 11.	reconstruct_species_phylogeny.py
    step_number = '11'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_species_phylogeny'
    script_path = os.path.join(consts.SRC_DIR, 'steps/reconstruct_species_phylogeny.py')
    phylogeny_path, phylogeny_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    start_time = time()
    if not os.path.exists(done_file_path):
        output_tree_path = os.path.join(phylogeny_path, 'final_species_tree.newick')
        if number_of_genomes < 3 or core_proteome_length == 0:
            logger.info('Skipping phylogeny reconstruction because there are less than 3 genomes or no core proteome')
            open(output_tree_path, 'w').close()
        else:
            logger.info('Reconstructing species phylogeny...')

            params = [aligned_core_proteome_file_path,
                      output_tree_path,
                      phylogeny_tmp_dir,
                      f'--cpu {consts.PHYLOGENY_NUM_OF_CORES}',
                      f'--bootstrap {"yes" if args.bootstrap else "no"}']
            if args.outgroup:
                with open(genomes_names_path, 'r') as genomes_names_fp:
                    strains_names = genomes_names_fp.read().split('\n')
                if args.outgroup in strains_names:
                    params += [f'--outgroup {args.outgroup}']
                else:
                    logger.info(f'Outgroup {args.outgroup} was specified but it is not one of the input species:\n'
                                f'{",".join(sorted(strains_names))}\nAn unrooted tree is going to be reconstructed')

            submit_mini_batch(logger, script_path, [params], phylogeny_tmp_dir, error_file_path,
                              args.queue_name, args.account_name, job_name='tree_reconstruction',
                              num_of_cpus=consts.PHYLOGENY_NUM_OF_CORES,
                              memory=consts.PHYLOGENY_REQUIRED_MEMORY_GB,
                              command_to_run_before_script='export QT_QPA_PLATFORM=offscreen')  # Needed to avoid an error in drawing the tree. Taken from: https://github.com/NVlabs/instant-ngp/discussions/300

            # wait for the phylogenetic tree here
            wait_for_results(logger, times_logger, step_name, phylogeny_tmp_dir,
                             num_of_expected_results=1, error_file_path=error_file_path,
                             start=start_time)

        add_results_to_final_dir(logger, phylogeny_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def step_12_codon_bias(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                       done_files_dir, orfs_dir, orthologs_dna_sequences_dir_path):
    # 12.  codon_bias.py
    # Input: ORF dir and OG dir
    # Output: W_vector for each genome, CAI for each OG
    step_number = '12'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_codon_bias'
    script_path = os.path.join(consts.SRC_DIR, 'steps/codon_bias.py')
    codon_bias_output_dir_path, codon_bias_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    cai_table_path = os.path.join(codon_bias_output_dir_path, 'CAI_table.csv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Analyzing codon bias...')
        params = [
            orfs_dir,
            orthologs_dna_sequences_dir_path,
            codon_bias_output_dir_path,
            cai_table_path,
            codon_bias_tmp_dir,
            consts.CODON_BIAS_NUM_OF_CORES
        ]
        submit_mini_batch(logger, script_path, [params], codon_bias_tmp_dir, error_file_path, args.queue_name, args.account_name, job_name='codon_bias',
                          num_of_cpus=consts.CODON_BIAS_NUM_OF_CORES)

        wait_for_results(logger, times_logger, step_name, codon_bias_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        add_results_to_final_dir(logger, codon_bias_output_dir_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return cai_table_path


def step_13_kegg_annotation(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                            done_files_dir, orthologs_aa_sequences_dir_path, final_orthologs_table_file_path):
    # 13.  kegg_annotation.py
    # Input: OG aa dir
    step_number = '13'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_kegg'
    script_path = os.path.join(consts.SRC_DIR, 'steps/kegg_annotation.py')
    output_dir_path, tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    kegg_table_path = os.path.join(output_dir_path, 'og_kegg.csv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Annotation with KEGG Orthology...')
        params = [
            orthologs_aa_sequences_dir_path,
            final_orthologs_table_file_path,
            output_dir_path,
            kegg_table_path,
            consts.KEGG_NUM_OF_CORES,
            '--optimize'
        ]
        submit_mini_batch(logger, script_path, [params], tmp_dir, error_file_path, args.queue_name, args.account_name,
                          job_name='kegg', num_of_cpus=consts.KEGG_NUM_OF_CORES, memory=consts.KEGG_REQUIRED_MEMORY_GB)

        wait_for_results(logger, times_logger, step_name, tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return kegg_table_path


def step_14_orthogroups_annotations(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                                    done_files_dir, final_orthologs_table_file_path, kegg_annotations, codon_bias_annotations):
    # 14.  add annotations to otrhogroups table
    step_number = '14'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_annotations'
    output_dir_path, tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Adding annotations to orthogroups...')
        final_orthologs_df = pd.read_csv(final_orthologs_table_file_path)

        if os.path.exists(kegg_annotations):
            kegg_table_df = pd.read_csv(kegg_annotations)[['OG_name', 'knum', 'knum_description']]
            final_orthologs_df = pd.merge(kegg_table_df, final_orthologs_df, on='OG_name')

        if os.path.exists(codon_bias_annotations):
            cai_df = pd.read_csv(codon_bias_annotations)[['OG_name', 'CAI_mean']]
            final_orthologs_df = pd.merge(cai_df, final_orthologs_df, on='OG_name')

        final_orthologs_table_annotated_path = os.path.join(output_dir_path, 'orthogroups_annotated.csv')
        final_orthologs_df.to_csv(final_orthologs_table_annotated_path, index=False)

        add_results_to_final_dir(logger, final_orthologs_table_annotated_path, final_output_dir,
                                 keep_in_source_dir=consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def run_main_pipeline(args, logger, times_logger, error_file_path, progressbar_file_path, output_html_path, output_dir,
                      tmp_dir, done_files_dir, data_path, genomes_names_path, final_output_dir):
    with open(genomes_names_path, 'r') as genomes_names_fp:
        genomes_names = genomes_names_fp.read().split('\n')
    number_of_genomes = len(genomes_names)

    if args.filter_out_plasmids or args.inputs_fasta_type == 'orfs':
        filtered_inputs_dir = step_0_fix_input_files(args, logger, times_logger, error_file_path, output_dir,
                                                     tmp_dir, done_files_dir, data_path)
        data_path = filtered_inputs_dir
        if args.filter_out_plasmids:
            update_progressbar(progressbar_file_path, 'Filter out plasmids')
        edit_progress(output_html_path, progress=5)

    if args.step_to_complete == '0':
        logger.info("Step 0 completed.")
        return

    if not args.only_calc_ogs:
        step_1_calculate_ani(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                             done_files_dir, data_path)
        update_progressbar(progressbar_file_path, 'Calculate ANI (Average Nucleotide Identity)')
        edit_progress(output_html_path, progress=10)

    if args.step_to_complete == '1':
        logger.info("Step 1 completed.")
        return

    orfs_dir, translated_orfs_dir, all_proteins_fasta_path = step_2_search_orfs(
        args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir, done_files_dir, data_path)
    update_progressbar(progressbar_file_path, 'Predict and translate ORFs')
    update_progressbar(progressbar_file_path, 'Translate ORFs')
    edit_progress(output_html_path, progress=15)

    if args.step_to_complete == '2':
        logger.info("Step 2 completed.")
        return

    if not args.only_calc_ogs:
        step_3_analyze_genome_completeness(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                           final_output_dir, done_files_dir, translated_orfs_dir)
        update_progressbar(progressbar_file_path, 'Calculate genomes completeness')
        edit_progress(output_html_path, progress=20)

    if args.step_to_complete == '3':
        logger.info("Step 3 completed.")
        return

    clusters_output_dir = step_4_cluster_proteomes(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                             done_files_dir, all_proteins_fasta_path)
    edit_progress(output_html_path, progress=35)

    if args.step_to_complete == '4':
        logger.info("Step 4 completed.")
        return

    orthologs_table_file_path = step_5_infer_orthogroups(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                             done_files_dir, clusters_output_dir, translated_orfs_dir, all_proteins_fasta_path, genomes_names_path)
    update_progressbar(progressbar_file_path, 'Infer orthogroups')
    edit_progress(output_html_path, progress=35)

    if args.step_to_complete == '5':
        logger.info("Step 5 completed.")
        return

    orphan_genes_dir = step_6_extract_orphan_genes(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                                   final_output_dir, done_files_dir, orfs_dir, orthologs_table_file_path)
    update_progressbar(progressbar_file_path, 'Find orphan genes')
    edit_progress(output_html_path, progress=35)

    if args.step_to_complete == '6':
        logger.info("Step 6 completed.")
        return

    final_orthologs_table_file_path = \
        step_7_orthologs_table_variations(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                          final_output_dir, done_files_dir, orthologs_table_file_path, orphan_genes_dir)
    edit_progress(output_html_path, progress=55)

    if args.step_to_complete == '7':
        logger.info("Step 7 completed.")
        return

    if args.only_calc_ogs:
        return

    ogs_dna_sequences_path, og_aa_sequences_path, ogs_aa_alignments_path, ogs_dna_alignments_path = \
        step_8_build_orthologous_groups_fastas(args, logger, times_logger, error_file_path,
                                               output_dir, tmp_dir, final_output_dir, done_files_dir, orfs_dir,
                                               final_orthologs_table_file_path)
    update_progressbar(progressbar_file_path, 'Prepare orthogroups fasta files')
    edit_progress(output_html_path, progress=60)

    if args.step_to_complete == '8':
        logger.info("Step 8 completed.")
        return

    aligned_core_proteome_file_path, core_proteome_length = step_9_extract_core_genome_and_core_proteome(
        args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir, done_files_dir,
        number_of_genomes, genomes_names_path, ogs_aa_alignments_path, ogs_dna_alignments_path)
    update_progressbar(progressbar_file_path, 'Infer core genome and proteome')
    edit_progress(output_html_path, progress=75)

    if args.step_to_complete == '9':
        logger.info("Step 9 completed.")
        return

    if args.core_minimal_percentage == 100:
        step_10_genome_numeric_representation(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                              final_output_dir, done_files_dir, orfs_dir, final_orthologs_table_file_path)
        update_progressbar(progressbar_file_path, 'Calculate genomes numeric representation')
        edit_progress(output_html_path, progress=90)

    if args.step_to_complete == '10':
        logger.info("Step 10 completed.")
        return

    step_11_phylogeny(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                      done_files_dir, aligned_core_proteome_file_path, genomes_names_path, number_of_genomes, core_proteome_length)
    update_progressbar(progressbar_file_path, 'Reconstruct species phylogeny')
    edit_progress(output_html_path, progress=95)

    if args.step_to_complete == '11':
        logger.info("Step 11 completed.")
        return

    cai_table_path = step_12_codon_bias(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                       done_files_dir, orfs_dir, ogs_dna_sequences_path)
    update_progressbar(progressbar_file_path, 'Analyze codon bias')
    edit_progress(output_html_path, progress=65)

    if args.step_to_complete == '12':
        logger.info("Step 12 completed.")
        return

    kegg_table_path = step_13_kegg_annotation(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                            done_files_dir, og_aa_sequences_path, final_orthologs_table_file_path)
    update_progressbar(progressbar_file_path, 'Annotate orthogroups with KEGG Orthology (KO) terms')
    edit_progress(output_html_path, progress=65)

    if args.step_to_complete == '13':
        logger.info("Step 13 completed.")
        return

    step_14_orthogroups_annotations(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                                    done_files_dir, final_orthologs_table_file_path, kegg_table_path, cai_table_path)
    edit_progress(output_html_path, progress=65)

    if args.step_to_complete == '14':
        logger.info("Step 14 completed.")
        return


def report_error_in_main_pipeline_to_admin(logger, e, meta_output_dir, error_file_path, run_number, output_html_path,
                                           meta_output_url):
    error_msg = f'{flask_interface_consts.WEBSERVER_NAME} failed :('
    if os.path.exists(error_file_path):
        with open(error_file_path) as error_f:
            error_txt = error_f.read()
            logger.error(f'error.txt file says:\n{error_txt}')
            error_msg = f'The job was failed due to the following reason:<br>{error_txt}'

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    logger.error(f'\n\n{"$" * 100}\n\n{error_msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\n'
                 f'e: {e}\n\n{traceback.format_exc()}\n\n{"$" * 100}')

    edit_failure_html(logger, output_html_path, run_number, error_msg)
    edit_progress(output_html_path, active=False)

    if consts.SEND_MAILS:
        notify_admin(meta_output_dir, meta_output_url, run_number)


def report_main_pipeline_result_to_user(args, logger, status, total_time, output_url, run_number):
    msg = f'M1CR0B1AL1Z3R pipeline {status}'
    if status == 'is done':
        msg += f' (Took {total_time}).\nResults can be found at {output_url}.\nPlease note that the ' \
               f'results will be kept in the server for three months.'
    else:
        msg += f'. For further information please visit: {output_url}'
    logger.info(msg)

    if consts.SEND_MAILS:
        logger.info(f'Sending a notification email to {args.email}')
        try:
            send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>', args.email,
                       subject=f'{flask_interface_consts.WEBSERVER_NAME} run number {run_number} {status}.', content=msg)
        except:
            logger.error(f'\nFailed sending notification to {args.email}\n')


# def run_pipeline_extensions(args, logger, times_logger, error_file_path, run_number, output_dir, tmp_dir,
#                             done_files_dir, data_path, orfs_dir,
#                             orfs_step_number, final_orthologs_table_file_path, phylogenetic_raw_tree_path,
#                             final_output_dir_name):
#     logger.info('\n\n\n')
#     logger.info('#' * 100)
#     logger.info('Executing additional modules for Sweeps analysis...')
#     # 21.	extract_promoters_and_orfs
#     # Input: (1) A path to a genome (2) Prodigal output file with ORFs coordinates
#     # Output: A fasta file containing promoters+genes (all ORFS of the given coordinates + k[=300] upstream bases)
#     # Can be parallelized on cluster
#     step_number = '21'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = f'{step_number}_extract_promoters_and_orfs'
#     script_path = os.path.join(consts.SRC_DIR, 'steps/extract_promoters_and_orfs.py')
#     pipeline_step_output_dir_21, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
#     done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
#     if not os.path.exists(done_file_path):
#         logger.info('Extracting promoters...')
#         all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
#         for genome_file in os.listdir(data_path):
#             genome_path = os.path.join(data_path, genome_file)
#             genome_file_prefix = os.path.splitext(genome_file)[0]
#             orfs_path = os.path.join(orfs_dir, f'{genome_file_prefix}.{orfs_step_number}_ORFs')
#             output_prefix = os.path.join(pipeline_step_output_dir_21, f'{genome_file_prefix}.promoter_and_orf')
#
#             single_cmd_params = [genome_path, orfs_path, output_prefix,
#                                  f'--promoters_length {args.promoters_length}']  # (optional)
#             all_cmds_params.append(single_cmd_params)
#
#         num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
#                                                    num_of_cmds_per_job=5,
#                                                    job_name_suffix='promoter_gene_alignment',
#                                                    queue_name=args.queue_name,
#                                                    required_modules_as_list=[consts.MAFFT])
#
#         wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
#                          num_of_batches, error_file_path, email=args.email)
#
#         write_to_file(logger, done_file_path, '.')
#     else:
#         logger.info(f'done file {done_file_path} already exists. Skipping step...')
#
#     # 22.	extract_orfs.py
#     # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
#     # Output: write the sequences of the orthologs group to the output file.
#     # Can be parallelized on cluster
#     step_number = '22'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = f'{step_number}_orthologs_groups_dna_sequences'
#     script_path = os.path.join(consts.SRC_DIR, 'steps/extract_orfs.py')
#     pipeline_step_output_dir_22, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
#     assert os.path.exists(
#         pipeline_step_output_dir_22), f'Failed to create output folder. {pipeline_step_output_dir_22} does not exist!'
#     done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
#     if not os.path.exists(done_file_path):
#         with open(final_orthologs_table_file_path) as f:
#             logger.info('Extracting orthologs groups sequences according to final orthologs table...')
#             all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
#             header_line = f.readline()
#             first_delimiter_index = header_line.index(consts.CSV_DELIMITER)
#             final_table_header = header_line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"
#             for line in f:
#                 first_delimiter_index = line.index(consts.CSV_DELIMITER)
#                 og_name = line[:first_delimiter_index]
#                 cluster_members = line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"
#
#                 number_of_sequences_in_cluster = sum(
#                     1 for member in cluster_members.split(consts.CSV_DELIMITER) if member)
#                 if number_of_sequences_in_cluster < args.minimal_number_of_sequences_allowed_for_sweeps_analysis:
#                     # ignore msa that contains less than the minimal allowed sequences for sweeps analysis
#                     logger.info(
#                         f'Orthologs group #{og_name} contains {number_of_sequences_in_cluster} members and thus '
#                         f'will not be included in the sweeps analysis (at least '
#                         f'{args.minimal_number_of_sequences_allowed_for_sweeps_analysis} '
#                         f'members are required)')
#                     continue
#
#                 logger.debug(f'Orthologs group #{og_name} contains {number_of_sequences_in_cluster} sequences and'
#                              f' thus will be included in the sweeps analysis')
#                 og_sequences_path = os.path.join(pipeline_step_output_dir_22, f'{og_name}_dna.fas')
#
#                 single_cmd_params = [pipeline_step_output_dir_21,
#                                      f'"{final_table_header}"',
#                                      # should be flanked by quotes because it might contain spaces...
#                                      f'"{cluster_members}"',
#                                      # should be flanked by quotes because it might contain spaces...
#                                      f'"{og_name}"',
#                                      # should be flanked by quotes because it might contain spaces...
#                                      og_sequences_path]
#                 all_cmds_params.append(single_cmd_params)
#
#         num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
#                                                    num_of_cmds_per_job=100,
#                                                    job_name_suffix='extract_sequences',
#                                                    queue_name=args.queue_name)
#
#         wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
#                          num_of_batches, error_file_path, email=args.email)
#
#         write_to_file(logger, done_file_path, '.')
#     else:
#         logger.info(f'done file {done_file_path} already exists. Skipping step...')
#
#     # 23.	align_orthologs_group.py
#     # Input: (1) A path to an unaligned sequences file (2) An output file path
#     # Output: aligned sequences
#     step_number = '23'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = f'{step_number}_aligned_dna_orthologs_groups_with_promoter'
#     script_path = os.path.join(consts.SRC_DIR, 'steps/align_orthologs_group.py')
#     pipeline_step_output_dir_23, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
#     done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
#     if not os.path.exists(done_file_path):
#         logger.info('Aligning orthologs groups...')
#         all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
#         for og_file in os.listdir(pipeline_step_output_dir_22):
#             og_path = os.path.join(pipeline_step_output_dir_22, og_file)
#             alignment_path = os.path.join(pipeline_step_output_dir_23, og_file.replace('.fas', '_mafft.fas'))
#
#             single_cmd_params = [og_path, alignment_path, '--type nuc']
#             all_cmds_params.append(single_cmd_params)
#
#         num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
#                                                    num_of_cmds_per_job=100,
#                                                    job_name_suffix='align_promoter_and_gene',
#                                                    queue_name=args.queue_name,
#                                                    required_modules_as_list=[consts.MAFFT])
#
#         wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
#                          num_of_batches, error_file_path, email=args.email)
#
#         write_to_file(logger, done_file_path, '.')
#     else:
#         logger.info(f'done file {done_file_path} already exists. Skipping step...')
#
#     # 24.	adjust_tree_to_msa.py
#     # Input: (1) A path to an MSA file to which the tree should be adjusted
#     #        (2) A path to a species tree that contains (at least) all the species in the input MSA
#     #        (3) A path to a folder in which a txt file (with the same name as the msa_file) will be created. All the species that do not appear in the msa (and thus will be removed) will be written to the file that was created',
#     #        (4) A path to a file in which the pruned tree will be written
#     # Output: an adjusted tree per msa at (4)
#     step_number = '24'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = f'{step_number}_prunned_trees'
#     script_path = f'/bioseq/sincopa/adjust_tree_to_msa.py'
#     pipeline_step_output_dir_24, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
#     done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
#     if not os.path.exists(done_file_path):
#         logger.info('Removing bootstrap values from tree...')
#         phylogenetic_raw_tree_path = phylogenetic_raw_tree_path.replace('outputs',
#                                                                         final_output_dir_name)  # the tree was moved to the final dir..
#         phylogenetic_raw_tree_path_without_bootstrap_values = phylogenetic_raw_tree_path.replace('.txt',
#                                                                                                  '_no_bootstrap.txt')
#         remove_bootstrap_values(phylogenetic_raw_tree_path, phylogenetic_raw_tree_path_without_bootstrap_values)
#
#         logger.info('Pruning trees...')
#         all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
#         for msa_file in os.listdir(pipeline_step_output_dir_23):
#             msa_path = os.path.join(pipeline_step_output_dir_23, msa_file)
#             pruned_tree_path = os.path.join(pipeline_step_output_dir_24, msa_file.replace('.fas', '.tree'))
#
#             single_cmd_params = [msa_path, phylogenetic_raw_tree_path_without_bootstrap_values,
#                                  pipeline_step_tmp_dir, pruned_tree_path]
#             all_cmds_params.append(single_cmd_params)
#
#         num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
#                                                    num_of_cmds_per_job=100,
#                                                    job_name_suffix='adjust_tree', queue_name=args.queue_name)
#
#         wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
#                          num_of_batches, error_file_path, email=args.email)
#
#         write_to_file(logger, done_file_path, '.')
#     else:
#         logger.info(f'done file {done_file_path} already exists. Skipping step...')
#
#     # 25.	fix_msa.py
#     # Input: (1) A path to an MSA file that (might) contain too few genomes and/or
#     #                                                       ambiguous characters (i.e., other than a, c, g, t, -)
#     #        (2) A path to a fixed msa
#     # Output: a fixed msa at (2)
#     step_number = '25'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = f'{step_number}_fixed_dna_msa'
#     pipeline_step_output_dir_25, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
#     done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
#     if not os.path.exists(done_file_path):
#         script_path = f'/bioseq/sincopa/fix_msa.py'
#         all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
#         for msa_file in os.listdir(pipeline_step_output_dir_23):
#             # preparing parameters
#             msa_path = os.path.join(pipeline_step_output_dir_23, msa_file)
#             fixed_msa_output_path = os.path.join(pipeline_step_output_dir_25, msa_file)
#
#             single_cmd_params = [msa_path, fixed_msa_output_path]
#             all_cmds_params.append(single_cmd_params)
#
#         num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
#                                                    num_of_cmds_per_job=250,
#                                                    job_name_suffix='fix_msa', queue_name=args.queue_name)
#
#         wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
#                          num_of_batches, error_file_path, email=args.email)
#
#         write_to_file(logger, done_file_path, '.')
#     else:
#         logger.info(f'done file {done_file_path} already exists. Skipping step...')
#
#     # 26.	compute_homoplasy.py
#     # Input: (1) A path to an MSA file without ambiguous characters (only 'a', 'c', 'g', 't', and '-' are allowed)
#     #        (2) A path to a corresponding background tree
#     #        (3) A path to a control file for the homoplasy calculation cpp script
#     #        (4) A path to a file in which the homoplasy will be written
#     # Output: homoplasy of the input msa at (2)
#     step_number = '26'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = f'{step_number}_homoplasy'
#     pipeline_step_output_dir_26, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
#     done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
#     if not os.path.exists(done_file_path):
#         script_path = f'/bioseq/sincopa/compute_homoplasy.py'
#         all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
#         for msa_file in os.listdir(pipeline_step_output_dir_23):
#             # preparing parameters
#             fixed_msa_path = os.path.join(pipeline_step_output_dir_25, msa_file)
#             tree_path = os.path.join(pipeline_step_output_dir_24, msa_file.replace('fas', 'tree'))
#             control_path = os.path.join(pipeline_step_tmp_dir, msa_file.replace('fas', 'control'))
#             output_path = os.path.join(pipeline_step_output_dir_26, msa_file.replace('fas', 'homoplasy'))
#
#             single_cmd_params = [fixed_msa_path, tree_path, control_path, output_path]
#             all_cmds_params.append(single_cmd_params)
#
#         num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
#                                                    job_name_suffix='homoplasy', queue_name=args.queue_name,
#                                                    num_of_cmds_per_job=50,
#                                                    required_modules_as_list=[consts.GCC])
#
#         wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
#                          num_of_batches, error_file_path, email=args.email)
#
#         write_to_file(logger, done_file_path, '.')
#     else:
#         logger.info(f'done file {done_file_path} already exists. Skipping step...')
#
#     # 27.	compute_sweeps_score.py
#     # Input: (1) A path to an MSA file to be analyzed
#     #        (2) A path to a homoplasy
#     #        (3) A path to a file in which the sweep scores will be written
#     #        (4) A path to a file in which the sweep scores plot will be saved
#     # Output: sweep scores and sweep scores plot at (3) and (4), respectively
#     step_number = '27'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = f'{step_number}_sweeps_scores_computation'
#     pipeline_step_output_dir_27, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
#     done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
#     window_size = 50
#     if not os.path.exists(done_file_path):
#         script_path = f'/bioseq/sincopa/compute_sweeps_score.py'
#         all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
#         for fixed_msa_file in os.listdir(pipeline_step_output_dir_25):
#             # preparing parameters
#             fixed_msa_path = os.path.join(pipeline_step_output_dir_25, fixed_msa_file)
#             homoplasy_path = os.path.join(pipeline_step_output_dir_26, fixed_msa_file.replace('fas', 'homoplasy'))
#             output_scores_path = os.path.join(pipeline_step_output_dir_27, fixed_msa_file.replace('fas', 'txt'))
#             output_meta_path = os.path.join(pipeline_step_output_dir_27, fixed_msa_file.replace('fas', 'csv'))
#             output_plot_path = os.path.join(pipeline_step_output_dir_27, fixed_msa_file.replace('fas', 'png'))
#
#             single_cmd_params = [fixed_msa_path, homoplasy_path, output_scores_path,
#                                  output_meta_path, output_plot_path, window_size]
#             all_cmds_params.append(single_cmd_params)
#
#         num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
#                                                    job_name_suffix='sweeps', queue_name=args.queue_name,
#                                                    num_of_cmds_per_job=100)
#
#         wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
#                          num_of_batches, error_file_path, email=args.email)
#
#         write_to_file(logger, done_file_path, '.')
#     else:
#         logger.info(f'done file {done_file_path} already exists. Skipping step...')
#
#     # 28. analyse sweeps
#     # Input: a path to a folder with sweeps score files
#     # Output: concatenated csv file with all the analysis metadata
#     # CANNOT be parallelized on cluster
#     step_number = '28'
#     logger.info(f'Step {step_number}: {"_" * 100}')
#     step_name = 'sweeps_analysis'
#     sweeps_analysis_dir = os.path.join(output_dir, f'{step_number}_{step_name}')
#     os.makedirs(sweeps_analysis_dir, exist_ok=True)
#     sweeps_summary_path = os.path.join(sweeps_analysis_dir, 'sweeps_summary.csv')
#     sorted_sweeps_summary_path = os.path.join(sweeps_analysis_dir, 'sorted_sweeps_summary.csv')
#
#     if not os.path.exists(sorted_sweeps_summary_path):
#         logger.info('Concatenating sweeps analysis metadata results...')
#
#         # write summary header
#         header = 'msa_name,max_score,number_of_sequences,centrality,msa_length,window_size,index_of_max,' \
#                  'mean_score,median_score,max_mean_division,max_median_division,relative_location_of_peak,' \
#                  'apd,pi,above95,above75,above50,above25,above05'
#         with open(sweeps_summary_path, 'w') as f:
#             f.write(f'{header}\n')
#
#         # add content
#         execute(logger, f'cat {pipeline_step_output_dir_27}/*csv >> {sweeps_summary_path}', process_is_string=True)
#
#         df = pd.read_csv(sweeps_summary_path)
#
#         df.sort_values(by=['max_score', 'number_of_sequences', 'centrality', 'msa_length'], ascending=False,
#                        inplace=True)
#         df.to_csv(sorted_sweeps_summary_path, index=False)
#
#         run_number2species = {'159716760987322741064374127401': 'C. trachomatis, r/m=0.66',
#                               '159717709089300979385014858842': 'H. pylori, r/m=21.11',
#                               '159717758554609942254414562338': 'S. pneumoniae, r/m=5.17',
#                               '159717782526726941828676365733': 'N. meningitidis, r/m=14.62',
#                               '159717817742964662930817004125': 'E. coli r/m=0.38',
#                               '159717825610100634511315581997': 'M. tuberculosis, r/m=0'}
#
#         for column in df.columns:
#             if column == 'msa_name' or 'above' in column:
#                 logger.info(f'Skipping {column} column...')
#                 continue
#
#             fig = plt.figure()
#             plt.title(f'{column}\n({run_number2species.get(run_number, "unrecognized species")})')
#             plt.hist(df[column], bins=50)  # TODO: density=True
#             fig.savefig(f'{sweeps_analysis_dir}/{column}.png', dpi=500, bbox_inches='tight')
#             plt.close()
#             fig = plt.figure()
#             plt.title(f'{column}\n({run_number2species.get(run_number, "unrecognized species")})')
#             df[column].plot.kde()
#             fig.savefig(f'{sweeps_analysis_dir}/{column}_smooth.png', dpi=500, bbox_inches='tight')
#             plt.close()
#             fig = plt.figure()
#             percentile = 10 / 100
#             plt.title(
#                 f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
#             df[column][df[column] >= df[column].quantile(1 - percentile)].plot.kde()
#             fig.savefig(f'{sweeps_analysis_dir}/{column}_smooth_top_{percentile * 100}_percentile.png', dpi=500,
#                         bbox_inches='tight')
#             plt.close()
#             fig = plt.figure()
#             percentile = 20 / 100
#             plt.title(
#                 f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
#             df[column][df[column] >= df[column].quantile(1 - percentile)].plot.kde()
#             fig.savefig(f'{sweeps_analysis_dir}/{column}_smooth_top_{percentile * 100}_percentile.png', dpi=500,
#                         bbox_inches='tight')
#             plt.close()
#             fig = plt.figure()
#             percentile = 5 / 1000
#             plt.title(
#                 f'{column} ({percentile} trimmed)\n({run_number2species.get(run_number, "unrecognized species")})')
#             plt.hist(np.sort(df[column])[int(len(df) * percentile): len(df) - int(len(df) * percentile)], bins=50)
#             fig.savefig(f'{sweeps_analysis_dir}/{column}_trimmed.png', dpi=500, bbox_inches='tight')
#             plt.close()
#             fig = plt.figure()
#             percentile = 1 / 100
#             plt.title(
#                 f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
#             plt.hist(np.sort(df[column])[len(df) - int(len(df) * percentile):], bins=30)
#             fig.savefig(f'{sweeps_analysis_dir}/{column}_top_{percentile * 100}_percentile.png', dpi=500,
#                         bbox_inches='tight')
#             plt.close()
#             fig = plt.figure()
#             percentile = 5 / 100
#             plt.title(
#                 f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
#             plt.hist(np.sort(df[column])[len(df) - int(len(df) * percentile):], bins=30)
#             fig.savefig(f'{sweeps_analysis_dir}/{column}_top_{percentile * 100}_percentile.png', dpi=500,
#                         bbox_inches='tight')
#             plt.close()
#
#         # No need to wait...
#     else:
#         logger.info(f'done file {sorted_sweeps_summary_path} already exists. Skipping step...')
#
#     # preventing folder's deletion (output dir is being deleted
#     # logger.info(f'Moving {pipeline_step_output_dir_24} TO {os.path.join(meta_output_dir, step_name)}')
#     # try:
#     #     shutil.move(pipeline_step_output_dir_23, meta_output_dir)
#     #     shutil.move(pipeline_step_output_dir_24, meta_output_dir)
#     # except FileExistsError:
#     #     pass


def main(args):
    start_time = time()

    logger, times_logger, meta_output_dir, error_file_path, progressbar_file_path, run_number, output_html_path, \
        output_url, meta_output_url, output_dir, tmp_dir, done_files_dir, steps_results_dir, data_path = \
        prepare_pipeline_framework(args)

    try:
        validate_arguments(args)
        initialize_progressbar(args, progressbar_file_path)

        done_file_path = os.path.join(done_files_dir, f'prepare_and_verify_inputs.txt')
        genomes_names_path = os.path.join(output_dir, 'genomes_names.txt')
        if not os.path.exists(done_file_path):
            prepare_and_verify_input_data(args, logger, meta_output_dir, error_file_path, data_path, genomes_names_path)
            write_to_file(logger, done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists. Skipping step...')

        update_progressbar(progressbar_file_path, 'Validate input files')

        final_output_dir_name = f'{flask_interface_consts.WEBSERVER_NAME}_{args.output_dir}'
        final_output_dir = os.path.join(meta_output_dir, final_output_dir_name)

        run_main_pipeline(args, logger, times_logger, error_file_path, progressbar_file_path,
                          output_html_path, steps_results_dir, tmp_dir, done_files_dir,
                          data_path, genomes_names_path, final_output_dir)

        if args.step_to_complete is None or args.step_to_complete == PIPELINE_STEPS[-1] or args.only_calc_ogs \
                or args.zip_results_in_partial_pipeline:
            logger.info('Zipping results folder...')
            shutil.make_archive(final_output_dir, 'zip', meta_output_dir, final_output_dir_name)
            update_progressbar(progressbar_file_path, 'Finalize results')

        logger.info('Editing results html...')
        edit_success_html(logger, output_html_path, meta_output_dir, final_output_dir_name, run_number)
        edit_progress(output_html_path, progress=100, active=False)

        # run_pipeline_extensions(args, logger, error_file_path, run_number, output_dir, tmp_dir, done_files_dir,
        #                         data_path,
        #                         orfs_dir, orfs_step_number, final_orthologs_table_file_path,
        #                         phylogenetic_raw_tree_path, final_output_dir_name)

        status = 'is done'

        # remove intermediate results
        if run_number.lower() != 'example' and not consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR:
            logger.info('Cleaning up intermediate results...')
            remove_path(logger, steps_results_dir)
            remove_path(logger, tmp_dir)
    except Exception as e:
        status = 'was failed'
        with open(error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
        report_error_in_main_pipeline_to_admin(logger, e, meta_output_dir, error_file_path, run_number,
                                               output_html_path,
                                               meta_output_url)

    total_time = measure_time(int(time() - start_time))
    times_logger.info(f'Total pipeline time: {total_time}')
    report_main_pipeline_result_to_user(args, logger, status, total_time, output_url, run_number)

    logger.info('Done.')


if __name__ == '__main__':
    print(f'sys.path is\n{sys.path}')
    arguments = get_arguments()

    main(arguments)
