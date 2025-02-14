import argparse
import logging
import math
import os
import shutil
import sys
from collections import defaultdict
from time import time
import traceback
import json
import subprocess
from datetime import timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from auxiliaries.file_writer import write_to_file
from auxiliaries.input_verifications import prepare_and_verify_input_data
from auxiliaries.pipeline_auxiliaries import (wait_for_results, prepare_directories, fail, submit_mini_batch,
                                              submit_batch, add_results_to_final_dir, remove_path, str_to_bool,
                                              send_email_in_pipeline_end)
from auxiliaries import consts
from flask import flask_interface_consts
from auxiliaries.logic_auxiliaries import (plot_genomes_histogram, update_progressbar, define_intervals)
from auxiliaries.cluster_mmseqs_hits_to_orthogroups import (cluster_mmseqs_hits_to_orthogroups,
                                                            unify_clusters_mmseqs_hits, run_mmseqs_and_extract_hits)
from flask.SharedConsts import USER_FILE_NAME_ZIP, USER_FILE_NAME_TAR, State

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
    parser.add_argument('--email', help='A notification will be sent once the pipeline is done to this email. Optional.')
    parser.add_argument('--job_name', help='Optional job name.')
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
    parser.add_argument('--node_name', help='The node name to submit jobs to', default=None)
    parser.add_argument('--step_to_complete', help='The final step to execute', default=None,
                        choices=[*PIPELINE_STEPS, None])
    parser.add_argument('--only_calc_ogs', help='Do only the necessary steps to calculate OGs', action='store_true')
    parser.add_argument('--bypass_number_of_genomes_limit', help='Bypass the limit on number of genomes',
                        action='store_true')
    parser.add_argument('--pre_cluster_orthogroups_inference', help='Optimize the orthogroups inference using heuristics',
                        action='store_true')
    parser.add_argument('--num_of_clusters_in_orthogroup_inference', help='Number of clusters in the optimization of '
                                                                          'the orthogroups inference. Relevant only if '
                                                                          'pre_cluster_orthogroups_inference is True',
                        default=10)
    parser.add_argument('--auto_pre_cluster', help='Decide number of clusters based on number of strains. When True, ignores pre_cluster_orthogroups_inference and num_of_clusters_in_orthogroup_inference', action='store_true')
    parser.add_argument('--unify_clusters_after_mmseqs', help='If True, unify the clusters after the extraction of hits.'
                                                              'Otherwise, unify the orthogroups at the end of inference.',
                        action='store_true')
    parser.add_argument('--run_optimized_mmseqs', help='Optimize the mmseqs run',
                        action='store_true')
    parser.add_argument('--use_parquet', help='When True, use parquet files when possible instead of csv',
                        action='store_true')
    parser.add_argument('--use_linux_to_parse_big_files', action='store_true')
    parser.add_argument('--mmseqs_use_dbs', action='store_true')
    parser.add_argument('--do_not_copy_outputs_to_final_results_dir', action='store_true')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

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
    formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)

    logger = logging.getLogger('main')
    main_file_handler = logging.FileHandler(os.path.join(output_dir, 'main_log.txt'), mode='a')
    main_file_handler.setFormatter(formatter)
    logger.addHandler(main_file_handler)
    logger.setLevel(level)

    times_logger = logging.getLogger('times')
    times_file_handler = logging.FileHandler(os.path.join(output_dir, 'times_log.txt'), mode='a')
    times_file_handler.setFormatter(formatter)
    times_logger.addHandler(times_file_handler)
    times_logger.setLevel(level)

    logger.info(args)
    logger.info(f'meta_output_dir is: {meta_output_dir}')
    logger.info(f'Created output_dir in: {output_dir}')

    final_output_dir_name = f'{flask_interface_consts.WEBSERVER_NAME}_{args.output_dir}'
    final_output_dir = os.path.join(meta_output_dir, final_output_dir_name)
    logger.info(f'Creating final_output_dir is: {final_output_dir}')
    os.makedirs(final_output_dir, exist_ok=True)

    error_file_path = os.path.join(final_output_dir, flask_interface_consts.ERROR_FILE_NAME)
    progressbar_file_path = os.path.join(meta_output_dir, flask_interface_consts.PROGRESSBAR_FILE_NAME)

    run_number = os.path.join(os.path.split(meta_output_dir)[1])
    logger.info(f'run_number is {run_number}')

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

    return (logger, times_logger, meta_output_dir, error_file_path, progressbar_file_path, run_number, output_dir,
            tmp_dir, done_files_dir, steps_results_dir, data_path, final_output_dir_name, final_output_dir)


def validate_arguments(args):
    if type(args.bootstrap) == str:
        args.bootstrap = str_to_bool(args.bootstrap)
    if type(args.filter_out_plasmids) == str:
        args.filter_out_plasmids = str_to_bool(args.filter_out_plasmids)
    if type(args.add_orphan_genes_to_ogs) == str:
        args.add_orphan_genes_to_ogs = str_to_bool(args.add_orphan_genes_to_ogs)
    if type(args.pre_cluster_orthogroups_inference) == str:
        args.pre_cluster_orthogroups_inference = str_to_bool(args.pre_cluster_orthogroups_inference)
    if type(args.unify_clusters_after_mmseqs) == str:
        args.unify_clusters_after_mmseqs = str_to_bool(args.unify_clusters_after_mmseqs)
    if type(args.use_parquet) == str:
        args.use_parquet = str_to_bool(args.use_parquet)
    if type(args.auto_pre_cluster) == str:
        args.auto_pre_cluster = str_to_bool(args.auto_pre_cluster)
    if type(args.mmseqs_use_dbs) == str:
        args.mmseqs_use_dbs = str_to_bool(args.mmseqs_use_dbs)
    if type(args.run_optimized_mmseqs) == str:
        args.run_optimized_mmseqs = str_to_bool(args.run_optimized_mmseqs)
    if type(args.use_linux_to_parse_big_files) == str:
        args.use_linux_to_parse_big_files = str_to_bool(args.use_linux_to_parse_big_files)
    if type(args.do_not_copy_outputs_to_final_results_dir) == str:
        args.do_not_copy_outputs_to_final_results_dir = str_to_bool(args.do_not_copy_outputs_to_final_results_dir)

    if args.pre_cluster_orthogroups_inference and args.unify_clusters_after_mmseqs and args.use_parquet:
        raise ValueError('If unify_clusters_after_mmseqs is True, we can not use parquet since we concat the csv hit files')

    if args.auto_pre_cluster and args.pre_cluster_orthogroups_inference:
        raise ValueError('auto_pre_cluster and pre_cluster_orthogroups_inference can not be both True')

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

    if not args.pre_cluster_orthogroups_inference:
        args.num_of_clusters_in_orthogroup_inference = 1

    args.num_of_clusters_in_orthogroup_inference = int(args.num_of_clusters_in_orthogroup_inference)
    if args.num_of_clusters_in_orthogroup_inference < 1:
        raise ValueError(f'num_of_clusters_in_orthogroup_inference argument {args.num_of_clusters_in_orthogroup_inference} has invalid value')


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


def step_1_fix_input_files(args, logger, times_logger, error_file_path, output_dir, tmp_dir, done_files_dir,
                           data_path):
    # 01.	drop_plasmids_and_fix_frames.py
    step_number = '01'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_fix_input_files'
    script_path = os.path.join(consts.SRC_DIR, 'steps/drop_plasmids_and_fix_frames.py')
    filtered_inputs_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Filtering plasmids out...')

        job_index_to_fasta_file_names = defaultdict(list)

        for i, fasta_file_name in enumerate(os.listdir(data_path)):
            job_index = i % consts.MAX_PARALLEL_JOBS
            job_index_to_fasta_file_names[job_index].append(fasta_file_name)

        jobs_inputs_dir = os.path.join(pipeline_step_tmp_dir, 'job_inputs')
        os.makedirs(jobs_inputs_dir, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_fasta_file_names in job_index_to_fasta_file_names.items():
            job_input_path = os.path.join(jobs_inputs_dir, f'{job_index}.txt')
            with open(job_input_path, 'w') as f:
                for fasta_file_name in job_fasta_file_names:
                    f.write(f'{os.path.join(data_path, fasta_file_name)}\n')

            single_cmd_params = [job_input_path, filtered_inputs_dir]

            if args.filter_out_plasmids:
                single_cmd_params.append('--drop_plasmids')
            if args.inputs_fasta_type == 'orfs':
                single_cmd_params.append('--fix_frames')

            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=1,
                                                   job_name_suffix='drop_plasmids',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return filtered_inputs_dir


def step_2_search_orfs(args, logger, times_logger, error_file_path,  output_dir, tmp_dir, final_output_dir,
                       done_files_dir, data_path):
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
    script_path = os.path.join(consts.SRC_DIR, 'steps/search_orfs.py')
    orfs_dir, orfs_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orfs_sequences_dir = os.path.join(orfs_dir, 'orfs_sequences')
    orfs_statistics_dir = os.path.join(orfs_dir, 'orfs_statistics')
    orfs_translated_dir = os.path.join(orfs_dir, 'orfs_translated')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting ORFs...')

        os.makedirs(orfs_sequences_dir, exist_ok=True)
        os.makedirs(orfs_statistics_dir, exist_ok=True)
        os.makedirs(orfs_translated_dir, exist_ok=True)

        job_index_to_fasta_file_names = defaultdict(list)

        for i, fasta_file_name in enumerate(os.listdir(data_path)):
            job_index = i % consts.MAX_PARALLEL_JOBS
            job_index_to_fasta_file_names[job_index].append(fasta_file_name)

        jobs_inputs_dir = os.path.join(orfs_tmp_dir, 'job_inputs')
        os.makedirs(jobs_inputs_dir, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_fasta_file_names in job_index_to_fasta_file_names.items():
            job_input_path = os.path.join(jobs_inputs_dir, f'{job_index}.txt')
            with open(job_input_path, 'w') as f:
                for fasta_file_name in job_fasta_file_names:
                    f.write(f'{os.path.join(data_path, fasta_file_name)}\n')

            single_cmd_params = [job_input_path, orfs_sequences_dir, orfs_statistics_dir, orfs_translated_dir,
                                 args.inputs_fasta_type]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, orfs_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=1,
                                                   job_name_suffix='search_orfs',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, orfs_tmp_dir, num_of_batches, error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orfs_sequences_dir, final_output_dir)
            add_results_to_final_dir(logger, orfs_translated_dir, final_output_dir)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 02_2.	plot_orfs_statistics
    step_number = '02_2'
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

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orfs_plots_path, final_output_dir)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 02_3.  Aggregate all protein fasta files to one file
    step_number = '02_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_concat_orfs'
    concat_orfs_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    all_orfs_fasta_path = os.path.join(concat_orfs_dir_path, 'all_orfs.fna')
    all_proteins_fasta_path = os.path.join(concat_orfs_dir_path, 'all_proteomes.faa')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Concatenating ORFs sequences...')

        cmd = f"cat {os.path.join(orfs_sequences_dir, '*')} > {all_orfs_fasta_path}"
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True)

        cmd = f"cat {os.path.join(orfs_translated_dir, '*')} > {all_proteins_fasta_path}"
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orfs_sequences_dir, orfs_translated_dir, all_orfs_fasta_path, all_proteins_fasta_path


def step_3_analyze_genome_completeness(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                       final_output_dir, done_files_dir, translated_orfs_dir_path):
    # 03.  assessing_genomes_completeness.py
    # Input: path to faa file
    # Output: calculate proteome completeness based on a dataset of profile-HMMs that represent core bacterial genes.
    step_number = '03'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_genomes_completeness'
    script_path = os.path.join(consts.SRC_DIR, 'steps/assessing_genome_completeness.py')
    genome_completeness_dir_path, genome_completeness_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Calculating genomes completeness...')

        job_index_to_fasta_file_names = defaultdict(list)

        for i, fasta_file_name in enumerate(os.listdir(translated_orfs_dir_path)):
            job_index = i % consts.MAX_PARALLEL_JOBS
            job_index_to_fasta_file_names[job_index].append(fasta_file_name)

        jobs_inputs_dir = os.path.join(genome_completeness_tmp_dir, 'job_inputs')
        os.makedirs(jobs_inputs_dir, exist_ok=True)
        genomes_output_dir_path = os.path.join(genome_completeness_dir_path, 'individual_proteomes_outputs')
        os.makedirs(genomes_output_dir_path, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_fasta_file_names in job_index_to_fasta_file_names.items():
            job_input_path = os.path.join(jobs_inputs_dir, f'{job_index}.txt')
            with open(job_input_path, 'w') as f:
                for fasta_file_name in job_fasta_file_names:
                    f.write(f'{os.path.join(translated_orfs_dir_path, fasta_file_name)}\n')

                single_cmd_params = [job_input_path, genomes_output_dir_path]
                all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, genome_completeness_tmp_dir,
                                                   error_file_path,
                                                   num_of_cmds_per_job=1,
                                                   job_name_suffix='genomes_completeness',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, genome_completeness_tmp_dir,
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

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, genome_completeness_dir_path, final_output_dir)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def step_4_cluster_proteomes(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                             done_files_dir, all_proteins_fasta_path, translated_orfs_dir):
    # 04_1.	cluster_proteomes.py
    step_number = '04_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_cluster_proteomes'
    script_path = os.path.join(consts.SRC_DIR, 'steps/cluster_proteomes.py')
    clusters_dir, clusters_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    clusters_file_path = os.path.join(clusters_dir, 'clusters.tsv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        cpus = consts.MMSEQS_CLUSTER_NUM_OF_CORES
        params = [all_proteins_fasta_path,
                  clusters_dir,
                  clusters_file_path,
                  consts.MMSEQS_CLUSTER_MIN_SEQ_ID,
                  consts.MMSEQS_CLUSTER_MIN_COVERAGE,
                  cpus,
                  f'--num_of_clusters_in_orthogroup_inference {args.num_of_clusters_in_orthogroup_inference}'
                  ]

        submit_mini_batch(logger, script_path, [params], clusters_tmp_dir, error_file_path, args.queue_name, args.account_name,
                          job_name='cluster_proteomes', num_of_cpus=cpus, memory=consts.MMSEQS_CLUSTER_REQUIRED_MEMORY_GB, node_name=args.node_name)
        wait_for_results(logger, times_logger, step_name, clusters_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 04_2.  filter_fasta_files_by_cluster.py
    step_number = '04_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_filter_fastas_by_cluster'
    script_path = os.path.join(consts.SRC_DIR, 'steps/filter_fasta_file.py')
    clusters_fasta_files, clusters_fasta_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
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
                params = [cluster_genes_path, all_proteins_fasta_path, os.path.join(cluster_fasta_files, 'all_proteomes.faa')]
                all_cmds_params.append(params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, clusters_fasta_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='filter_fastas_by_clusters',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, clusters_fasta_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return clusters_fasta_files


def step_5_infer_orthogroups_clustered(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                       done_files_dir, clusters_fasta_files, genomes_names_path):
    # infer_orthogroups.py
    step_number = '05_1' if args.unify_clusters_after_mmseqs else '05'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_infer_orthogroups'
    script_path = os.path.join(consts.SRC_DIR, 'steps/infer_orthogroups.py')
    infer_orthogroups_output_dir, infer_orthogroups_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    infer_orthogroups_done_dir = os.path.join(done_files_dir, step_name)
    os.makedirs(infer_orthogroups_done_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs

        for cluster_fastas_dir_name in os.listdir(clusters_fasta_files):
            cluster_fastas_dir_path = os.path.join(clusters_fasta_files, cluster_fastas_dir_name)
            cluster_id = cluster_fastas_dir_name.split('_')[1]

            cluster_output_dir = os.path.join(infer_orthogroups_output_dir, str(cluster_id))
            os.makedirs(cluster_output_dir, exist_ok=True)
            cluster_tmp_dir = os.path.join(infer_orthogroups_tmp_dir, str(cluster_id))
            os.makedirs(cluster_tmp_dir, exist_ok=True)
            cluster_done_dir = os.path.join(infer_orthogroups_done_dir, str(cluster_id))
            os.makedirs(cluster_done_dir, exist_ok=True)

            params = [step_number, cluster_output_dir, cluster_tmp_dir, cluster_done_dir, cluster_fastas_dir_path,
                      os.path.join(cluster_fastas_dir_path, 'all_proteomes.faa'), genomes_names_path, args.queue_name,
                      args.account_name, args.node_name, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff,
                      max(1, consts.MAX_PARALLEL_JOBS // args.num_of_clusters_in_orthogroup_inference)]

            if args.run_optimized_mmseqs:
                params.append('--run_optimized_mmseqs')
            if args.unify_clusters_after_mmseqs:
                params.append('--unify_clusters_after_mmseqs')
            if args.use_parquet:
                params.append('--use_parquet')
            if args.verbose:
                params.append('--verbose')
            if args.use_linux_to_parse_big_files:
                params.append('--use_linux_to_parse_big_files')
            if args.mmseqs_use_dbs:
                params.append('--mmseqs_use_dbs')
            all_cmds_params.append(params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, infer_orthogroups_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='infer_orthogroups',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   time_in_hours=consts.INFER_ORTHOGROUPS_JOB_TIME_LIMIT_HOURS,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, infer_orthogroups_tmp_dir,
                         num_of_batches, error_file_path, recursive_step=True)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    if args.unify_clusters_after_mmseqs:
        orthologs_output_dir, orthologs_scores_statistics_dir, paralogs_output_dir, paralogs_scores_statistics_dir = \
            unify_clusters_mmseqs_hits(logger, times_logger, output_dir, tmp_dir, done_files_dir, error_file_path,
                                       infer_orthogroups_output_dir, args.run_optimized_mmseqs, args.queue_name,
                                       args.account_name, args.node_name, '05', 2)
        orthogroups_file_path = \
            cluster_mmseqs_hits_to_orthogroups(logger, times_logger, error_file_path, output_dir, tmp_dir,
                                               done_files_dir, orthologs_output_dir, orthologs_scores_statistics_dir,
                                               paralogs_output_dir, paralogs_scores_statistics_dir,
                                               consts.MAX_PARALLEL_JOBS, '05', 4, args.account_name, args.queue_name, args.node_name,
                                               args.use_parquet, genomes_names_path)
    else:  # Aggregate OG tables of all clusters
        all_og_tables = []
        for cluster_dir_name in os.listdir(infer_orthogroups_output_dir):
            cluster_og_table = os.path.join(infer_orthogroups_output_dir, cluster_dir_name, '05_10_verified_table', 'orthogroups.csv')
            cluster_og_df = pd.read_csv(cluster_og_table)
            all_og_tables.append(cluster_og_df)
        unified_og_table = pd.concat(all_og_tables, axis=0, ignore_index=True)
        unified_og_table['OG_name'] = [f'OG_{i}' for i in range(len(unified_og_table.index))]
        orthogroups_file_path = os.path.join(infer_orthogroups_output_dir, 'orthogroups.csv')
        unified_og_table.to_csv(orthogroups_file_path, index=False)

    return orthogroups_file_path


def step_5_infer_orthogroups_non_clustered(args, logger, times_logger, error_file_path, output_dir, tmp_dir, done_files_dir,
                                           translated_orfs_dir, all_proteins_path, strains_names_path):
    orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir = \
        run_mmseqs_and_extract_hits(logger, times_logger, '05', error_file_path, output_dir, tmp_dir,
                                    done_files_dir, translated_orfs_dir, all_proteins_path, strains_names_path,
                                    args.queue_name, args.account_name, args.node_name, args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff,
                                    consts.MAX_PARALLEL_JOBS, args.run_optimized_mmseqs, args.use_parquet,
                                    args.use_linux_to_parse_big_files, args.mmseqs_use_dbs, args.verbose)

    orthogroups_file_path = \
        cluster_mmseqs_hits_to_orthogroups(logger, times_logger, error_file_path, output_dir, tmp_dir,
                                           done_files_dir, orthologs_output_dir, orthologs_scores_statistics_dir,
                                           paralogs_output_dir, paralogs_scores_statistics_dir,
                                           consts.MAX_PARALLEL_JOBS, '05', 4, args.account_name, args.queue_name, args.node_name,
                                           args.use_parquet, strains_names_path)

    return orthogroups_file_path


def step_5_full_orthogroups_inference(args, logger, times_logger, error_file_path, output_dir, tmp_dir, done_files_dir,
                                      translated_orfs_dir, all_proteins_fasta_path, genomes_names_path):
    if args.pre_cluster_orthogroups_inference:
        clusters_output_dir = step_4_cluster_proteomes(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                                       done_files_dir, all_proteins_fasta_path, translated_orfs_dir)

        if args.step_to_complete == '4':
            logger.info("Step 4 completed.")
            return

        orthogroups_file_path = step_5_infer_orthogroups_clustered(args, logger, times_logger, error_file_path,
                                                                   output_dir, tmp_dir, done_files_dir,
                                                                   clusters_output_dir, genomes_names_path)
    else:
        orthogroups_file_path = step_5_infer_orthogroups_non_clustered(args, logger, times_logger, error_file_path,
                                                                       output_dir,
                                                                       tmp_dir, done_files_dir, translated_orfs_dir,
                                                                       all_proteins_fasta_path, genomes_names_path)

    return orthogroups_file_path


def step_6_extract_orphan_genes(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                                done_files_dir, orfs_dir, orthologs_table_file_path):
    # 6.   extract_orphan_genes.py
    step_number = '06'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orphan_genes'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_orphan_genes.py')
    orphan_genes_dir, orphans_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orphan_genes_internal_dir = os.path.join(orphan_genes_dir, 'orphans_lists_per_genome')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting orphan genes...')
        os.makedirs(orphan_genes_internal_dir, exist_ok=True)

        job_index_to_fasta_file_names = defaultdict(list)

        for i, fasta_file_name in enumerate(os.listdir(orfs_dir)):
            job_index = i % consts.MAX_PARALLEL_JOBS
            job_index_to_fasta_file_names[job_index].append(fasta_file_name)

        jobs_inputs_dir = os.path.join(orphans_tmp_dir, 'job_inputs')
        os.makedirs(jobs_inputs_dir, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for job_index, job_fasta_file_names in job_index_to_fasta_file_names.items():
            job_input_path = os.path.join(jobs_inputs_dir, f'{job_index}.txt')
            with open(job_input_path, 'w') as f:
                for fasta_file_name in job_fasta_file_names:
                    f.write(f'{os.path.join(orfs_dir, fasta_file_name)}\n')

            single_cmd_params = [job_input_path, orthologs_table_file_path, orphan_genes_internal_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, orphans_tmp_dir, error_file_path,
                                                   num_of_cmds_per_job=max(1, len(all_cmds_params) // consts.MAX_PARALLEL_JOBS),
                                                   job_name_suffix='extract_orphans',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, orphans_tmp_dir,
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

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orphan_genes_dir, final_output_dir)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orphan_genes_internal_dir


def step_7_orthologs_table_variations(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                      final_output_dir, done_files_dir, orthologs_table_file_path, orphan_genes_dir):
    # 07_1.	orthologs_table_variations.py
    # Input: (1) path to the orthologs table (2) path to the orphan genes directory.
    # Output: build variations of the orthologs table.
    step_number = '07_1'
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
                          job_name='final_ortholog_groups', node_name=args.node_name)
        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, final_orthologs_table_dir_path, final_output_dir)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 07_2.	extract_groups_sizes_frequency
    step_number = '07_2'
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
        plt.title('Orthogroups sizes distribution', fontsize=20, loc='center', wrap=True)
        plt.xlabel('OG size (number of genomes)', fontsize=15)
        plt.ylabel('Count of OGs of each OG size', fontsize=15)
        plt.tight_layout()
        plt.savefig(os.path.join(group_sizes_path, 'groups_sizes.png'), dpi=600)
        plt.clf()

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, group_sizes_path, final_output_dir)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return final_orthologs_table_file_path


def step_8_build_orthologous_groups_fastas(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                           final_output_dir, done_files_dir, all_orfs_fasta_path,
                                           all_proteins_fasta_path, final_orthologs_table_file_path):
    # 08_1.	extract_orfs.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    step_number = '08_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_fasta'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_orfs.py')
    orthogroups_fasta_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)

    orthogroups_dna_dir_path = os.path.join(orthogroups_fasta_dir_path, 'orthogroups_dna')
    orthologs_aa_dir_path = os.path.join(orthogroups_fasta_dir_path, 'orthogroups_aa')
    orthogroups_aa_msa_dir_path = os.path.join(orthogroups_fasta_dir_path, 'orthogroups_aa_msa')
    orthogroups_induced_dna_msa_dir_path = os.path.join(orthogroups_fasta_dir_path, 'orthogroups_induced_dna_msa')

    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')

        os.makedirs(orthogroups_dna_dir_path, exist_ok=True)
        os.makedirs(orthologs_aa_dir_path, exist_ok=True)
        os.makedirs(orthogroups_aa_msa_dir_path, exist_ok=True)
        os.makedirs(orthogroups_induced_dna_msa_dir_path, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        with open(final_orthologs_table_file_path, 'r') as fp:
            number_of_ogs = sum(1 for _ in fp) - 1

        ogs_intervals = define_intervals(0, number_of_ogs - 1, consts.MAX_PARALLEL_JOBS)
        for og_number_start, og_number_end in ogs_intervals:
            single_cmd_params = [all_orfs_fasta_path,
                                 all_proteins_fasta_path,
                                 final_orthologs_table_file_path,
                                 og_number_start,
                                 og_number_end,
                                 orthogroups_dna_dir_path,
                                 orthologs_aa_dir_path,
                                 orthogroups_aa_msa_dir_path,
                                 orthogroups_induced_dna_msa_dir_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   error_file_path,
                                                   num_of_cmds_per_job=1,
                                                   job_name_suffix='orfs_extraction',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, num_of_batches, error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, orthogroups_dna_dir_path, final_output_dir)
            add_results_to_final_dir(logger, orthologs_aa_dir_path, final_output_dir)
            add_results_to_final_dir(logger, orthogroups_aa_msa_dir_path, final_output_dir)
            add_results_to_final_dir(logger, orthogroups_induced_dna_msa_dir_path, final_output_dir)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthogroups_dna_dir_path, orthologs_aa_dir_path, orthogroups_aa_msa_dir_path, orthogroups_induced_dna_msa_dir_path


def step_9_extract_core_genome_and_core_proteome(args, logger, times_logger, error_file_path,
                                                  output_dir, tmp_dir, final_output_dir, done_files_dir,
                                                  strains_names_path, aa_alignments_path, dna_alignments_path):
    # 09_1.	extract aligned_core_proteome.py
    step_number = '09_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    core_proteome_step_name = f'{step_number}_aligned_core_proteome'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_core_genome.py')
    aligned_core_proteome_path, aligned_core_proteome_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, core_proteome_step_name)
    core_proteome_done_file_path = os.path.join(done_files_dir, f'{core_proteome_step_name}.txt')

    if not os.path.exists(core_proteome_done_file_path):
        logger.info('Extracting aligned core proteome...')

        aligned_core_proteome_file_path = os.path.join(aligned_core_proteome_path, 'aligned_core_proteome.fas')
        core_proteome_length_file_path = os.path.join(aligned_core_proteome_path, 'core_length.txt')

        params = [aa_alignments_path, strains_names_path,
                  aligned_core_proteome_file_path,
                  core_proteome_length_file_path,
                  f'--core_minimal_percentage {args.core_minimal_percentage}']  # how many members induce a core group?
        submit_mini_batch(logger, script_path, [params], aligned_core_proteome_tmp_dir, error_file_path,
                          args.queue_name, args.account_name, job_name='core_proteome', node_name=args.node_name)

    else:
        logger.info(f'done file {core_proteome_done_file_path} already exists. Skipping step...')


    # 09_2.      extract_aligned_core_genome
    step_number = '09_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    core_genome_step_name = f'{step_number}_aligned_core_genome'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_core_genome.py')
    aligned_core_genome_path, aligned_core_genome_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, core_genome_step_name)
    core_genome_done_file_path = os.path.join(done_files_dir, f'{core_genome_step_name}.txt')
    if not os.path.exists(core_genome_done_file_path):
        logger.info('Extracting aligned core genome...')

        aligned_core_genome_file_path = os.path.join(aligned_core_genome_path, 'aligned_core_genome.fas')
        core_genome_length_file_path = os.path.join(aligned_core_genome_path, 'core_length.txt')

        params = [dna_alignments_path, strains_names_path,
                  aligned_core_genome_file_path,
                  core_genome_length_file_path,
                  f'--core_minimal_percentage {args.core_minimal_percentage}']  # how many members induce a core group?
        submit_mini_batch(logger, script_path, [params], aligned_core_genome_tmp_dir, error_file_path,
                          args.queue_name, args.account_name, job_name='core_genome', node_name=args.node_name)

    else:
        logger.info(f'done file {core_genome_done_file_path} already exists. Skipping step...')


    # 09_3.	extract aligned_core_proteome.py (for phylogeny reconstruction)
    step_number = '09_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    core_proteome_reduced_step_name = f'{step_number}_aligned_core_proteome_reduced'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_core_genome.py')
    aligned_core_proteome_reduced_path, aligned_core_proteome_reduced_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, core_proteome_reduced_step_name)
    aligned_core_proteome_reduced_file_path = os.path.join(aligned_core_proteome_reduced_path, 'aligned_core_proteome.fas')
    core_proteome_reduced_length_file_path = os.path.join(aligned_core_proteome_reduced_path, 'core_length.txt')
    core_proteome_reduced_done_file_path = os.path.join(done_files_dir, f'{core_proteome_reduced_step_name}.txt')
    if not os.path.exists(core_proteome_reduced_done_file_path):
        logger.info('Extracting aligned core proteome for phylogeny reconstruction...')

        params = [aa_alignments_path, strains_names_path,
                  aligned_core_proteome_reduced_file_path,
                  core_proteome_reduced_length_file_path,
                  f'--core_minimal_percentage {args.core_minimal_percentage}',  # how many members induce a core group?
                  f'--max_number_of_ogs {consts.MAX_NUMBER_OF_CORE_OGS_FOR_PHYLOGENY}']
        submit_mini_batch(logger, script_path, [params], aligned_core_proteome_reduced_tmp_dir, error_file_path,
                          args.queue_name, args.account_name, job_name='core_proteome', node_name=args.node_name)

    else:
        logger.info(f'done file {core_proteome_reduced_done_file_path} already exists. Skipping step...')

    # Wait for all results here (since the steps aren't dependent on each other)
    if not os.path.exists(core_proteome_done_file_path):
        wait_for_results(logger, times_logger, core_proteome_step_name, aligned_core_proteome_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, aligned_core_proteome_path, final_output_dir)
        write_to_file(logger, core_proteome_done_file_path, '.')

    if not os.path.exists(core_genome_done_file_path):
        wait_for_results(logger, times_logger, core_genome_step_name, aligned_core_genome_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, aligned_core_genome_path, final_output_dir)
        write_to_file(logger, core_genome_done_file_path, '.')

    if not os.path.exists(core_proteome_reduced_done_file_path):
        wait_for_results(logger, times_logger, core_proteome_reduced_step_name, aligned_core_proteome_reduced_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)
        write_to_file(logger, core_proteome_reduced_done_file_path, '.')

    with open(core_proteome_reduced_length_file_path, 'r') as fp:
        core_proteome_reduced_length = int(fp.read().strip())

    return aligned_core_proteome_reduced_file_path, core_proteome_reduced_length


def step_10_genome_numeric_representation(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                          final_output_dir, done_files_dir, orfs_dir, final_orthologs_table_file_path):
    # 10.	genome_numeric_representation.py
    step_number = '10'
    logger.info(f'Step {step_number}: {"_" * 100}')
    numeric_step_name = f'{step_number}_genome_numeric_representation'
    script_path = os.path.join(consts.SRC_DIR, 'steps/genome_numeric_representation.py')
    numeric_representation_output_dir, numeric_representation_tmp_dir = prepare_directories(logger, output_dir,
                                                                                            tmp_dir, numeric_step_name)
    core_genome_numeric_representation_file_path = os.path.join(numeric_representation_output_dir,
                                                                'core_genome_numeric_representation.txt')
    numeric_done_file_path = os.path.join(done_files_dir, f'{numeric_step_name}.txt')
    if not os.path.exists(numeric_done_file_path):
        params = [final_orthologs_table_file_path,
                  orfs_dir,
                  core_genome_numeric_representation_file_path,
                  numeric_representation_tmp_dir
                  ]
        submit_mini_batch(logger, script_path, [params], numeric_representation_tmp_dir, error_file_path,
                          args.queue_name, args.account_name, job_name='numeric_representation', node_name=args.node_name)

        wait_for_results(logger, times_logger, numeric_step_name, numeric_representation_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, numeric_representation_output_dir, final_output_dir)
        write_to_file(logger, numeric_done_file_path, '.')
    else:
        logger.info(f'done file {numeric_done_file_path} already exists. Skipping step...')


def step_11_phylogeny(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                      done_files_dir, aligned_core_proteome_file_path, genomes_names_path, number_of_genomes,
                      core_proteome_length, inputs_data_path):
    # 11_1.   ani.py
    step_number = '11_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    ani_step_name = f'{step_number}_ani'
    script_path = os.path.join(consts.SRC_DIR, 'steps/ani.py')
    ani_output_dir, ani_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, ani_step_name)
    ani_done_file_path = os.path.join(done_files_dir, f'{ani_step_name}.txt')
    if not os.path.exists(ani_done_file_path):
        logger.info('Calculating ANI values...')
        ani_tmp_files = os.path.join(ani_tmp_dir, 'temp_results')
        os.makedirs(ani_tmp_files, exist_ok=True)

        genomes_paths = [os.path.join(inputs_data_path, genome_file_name) for genome_file_name in os.listdir(inputs_data_path)]
        genomes_list_path = os.path.join(ani_tmp_files, 'genomes_list.txt')
        with open(genomes_list_path, 'w') as genomes_list_file:
            genomes_list_file.write('\n'.join(genomes_paths))

        single_cmd_params = [genomes_list_path, ani_output_dir, f'--cpus {consts.ANI_NUM_OF_CORES}']

        submit_mini_batch(logger, script_path, [single_cmd_params], ani_tmp_dir, error_file_path, args.queue_name, args.account_name,
                          job_name='ANI', num_of_cpus=consts.ANI_NUM_OF_CORES, memory=consts.ANI_REQUIRED_MEMORY_GB, node_name=args.node_name)

    else:
        logger.info(f'done file {ani_done_file_path} already exists. Skipping step...')

    # 11_2.	reconstruct_species_phylogeny.py
    step_number = '11_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    phylogeny_step_name = f'{step_number}_species_phylogeny'
    script_path = os.path.join(consts.SRC_DIR, 'steps/reconstruct_species_phylogeny.py')
    phylogeny_path, phylogeny_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, phylogeny_step_name)
    phylogeny_done_file_path = os.path.join(done_files_dir, f'{phylogeny_step_name}.txt')
    start_time = time()
    if not os.path.exists(phylogeny_done_file_path):
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

            # Needed to avoid an error in drawing the tree (since the default xdg_runtime_dir is sometimes not writable)
            xdg_runtime_dir = os.path.join(phylogeny_tmp_dir, 'xdg_runtime_dir')
            os.makedirs(xdg_runtime_dir, exist_ok=True)
            os.chmod(xdg_runtime_dir, 0o700)

            submit_mini_batch(logger, script_path, [params], phylogeny_tmp_dir, error_file_path,
                              args.queue_name, args.account_name, job_name='tree_reconstruction',
                              num_of_cpus=consts.PHYLOGENY_NUM_OF_CORES,
                              memory=consts.PHYLOGENY_REQUIRED_MEMORY_GB,
                              command_to_run_before_script=f'export QT_QPA_PLATFORM=offscreen\nexport XDG_RUNTIME_DIR={xdg_runtime_dir}', # Needed to avoid an error in drawing the tree. Taken from: https://github.com/NVlabs/instant-ngp/discussions/300
                              node_name=args.node_name,
                              time_in_hours=consts.PHYLOGENY_JOB_TIME_LIMIT_HOURS)

            # wait for the phylogenetic tree here
            wait_for_results(logger, times_logger, phylogeny_step_name, phylogeny_tmp_dir,
                             num_of_expected_results=1, error_file_path=error_file_path,
                             start=start_time)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, phylogeny_path, final_output_dir)
        write_to_file(logger, phylogeny_done_file_path, '.')
    else:
        logger.info(f'done file {phylogeny_done_file_path} already exists. Skipping step...')

    # Wait for ANI here, such that ANI and phylogeny will run in parallel
    if not os.path.exists(ani_done_file_path):
        wait_for_results(logger, times_logger, ani_step_name, ani_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, ani_output_dir, final_output_dir)
        write_to_file(logger, ani_done_file_path, '.')


def step_12_orthogroups_annotations(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                                    done_files_dir, orfs_dir, orthologs_dna_sequences_dir_path, orthologs_aa_sequences_dir_path,
                                    final_orthologs_table_file_path):
    # 12_1.  codon_bias.py
    # Input: ORF dir and OG dir
    # Output: W_vector for each genome, CAI for each OG
    step_number = '12_1'
    logger.info(f'Step {step_number}: {"_" * 100}')
    codon_bias_step_name = f'{step_number}_codon_bias'
    script_path = os.path.join(consts.SRC_DIR, 'steps/codon_bias.py')
    codon_bias_output_dir_path, codon_bias_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, codon_bias_step_name)
    cai_table_path = os.path.join(codon_bias_output_dir_path, 'CAI_table.csv')
    codon_bias_done_file_path = os.path.join(done_files_dir, f'{codon_bias_step_name}.txt')
    if not os.path.exists(codon_bias_done_file_path):
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
                          num_of_cpus=consts.CODON_BIAS_NUM_OF_CORES, node_name=args.node_name)

    else:
        logger.info(f'done file {codon_bias_done_file_path} already exists. Skipping step...')

    # 12_2.  kegg_annotation.py
    # Input: OG aa dir
    step_number = '12_2'
    logger.info(f'Step {step_number}: {"_" * 100}')
    kegg_step_name = f'{step_number}_kegg'
    script_path = os.path.join(consts.SRC_DIR, 'steps/kegg_annotation.py')
    kegg_output_dir_path, kegg_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, kegg_step_name)
    kegg_table_path = os.path.join(kegg_output_dir_path, 'og_kegg.csv')
    kegg_done_file_path = os.path.join(done_files_dir, f'{kegg_step_name}.txt')
    if not os.path.exists(kegg_done_file_path):
        logger.info('Annotation with KEGG Orthology...')
        params = [
            orthologs_aa_sequences_dir_path,
            final_orthologs_table_file_path,
            kegg_output_dir_path,
            kegg_table_path,
            consts.KEGG_NUM_OF_CORES,
            '--optimize'
        ]
        submit_mini_batch(logger, script_path, [params], kegg_tmp_dir, error_file_path, args.queue_name, args.account_name,
                          job_name='kegg', num_of_cpus=consts.KEGG_NUM_OF_CORES, memory=consts.KEGG_REQUIRED_MEMORY_GB, node_name=args.node_name)

    else:
        logger.info(f'done file {kegg_done_file_path} already exists. Skipping step...')

    # Wait for results here (since the steps aren't dependent on each other)
    if not os.path.exists(codon_bias_done_file_path):
        wait_for_results(logger, times_logger, codon_bias_step_name, codon_bias_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, codon_bias_output_dir_path, final_output_dir)
        write_to_file(logger, codon_bias_done_file_path, '.')

    if not os.path.exists(kegg_done_file_path):
        wait_for_results(logger, times_logger, kegg_step_name, kegg_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path)
        write_to_file(logger, kegg_done_file_path, '.')

    # 12_3.  add annotations to otrhogroups table
    step_number = '12_3'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_annotations'
    output_dir_path, tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Adding annotations to orthogroups...')
        final_orthologs_df = pd.read_csv(final_orthologs_table_file_path)

        if os.path.exists(kegg_table_path):
            kegg_table_df = pd.read_csv(kegg_table_path)[['OG_name', 'knum', 'knum_description']]
            final_orthologs_df = pd.merge(kegg_table_df, final_orthologs_df, on='OG_name')

        if os.path.exists(cai_table_path):
            cai_df = pd.read_csv(cai_table_path)[['OG_name', 'CAI_mean']]
            final_orthologs_df = pd.merge(cai_df, final_orthologs_df, on='OG_name')

        final_orthologs_table_annotated_path = os.path.join(output_dir_path, 'orthogroups_annotated.csv')
        final_orthologs_df.to_csv(final_orthologs_table_annotated_path, index=False)
        logger.info(f'Final orthologs table with annotations saved to {final_orthologs_table_annotated_path}')

        if not args.do_not_copy_outputs_to_final_results_dir:
            add_results_to_final_dir(logger, final_orthologs_table_annotated_path, final_output_dir)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def run_main_pipeline(args, logger, times_logger, error_file_path, progressbar_file_path, output_dir,
                      tmp_dir, done_files_dir, data_path, genomes_names_path, final_output_dir):
    with open(genomes_names_path, 'r') as genomes_names_fp:
        genomes_names = genomes_names_fp.read().split('\n')
    number_of_genomes = len(genomes_names)

    if args.auto_pre_cluster and number_of_genomes > 50:
        args.pre_cluster_orthogroups_inference = True
        args.num_of_clusters_in_orthogroup_inference = int(math.sqrt(number_of_genomes))

    if args.filter_out_plasmids or args.inputs_fasta_type == 'orfs':
        filtered_inputs_dir = step_1_fix_input_files(args, logger, times_logger, error_file_path, output_dir,
                                                     tmp_dir, done_files_dir, data_path)
        data_path = filtered_inputs_dir
        if args.filter_out_plasmids:
            update_progressbar(progressbar_file_path, 'Filter out plasmids')

    if args.step_to_complete == '1':
        logger.info("Step 1 completed.")
        return

    orfs_dir, translated_orfs_dir, all_orfs_fasta_path, all_proteins_fasta_path = step_2_search_orfs(
        args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir, done_files_dir, data_path)
    update_progressbar(progressbar_file_path, 'Predict and translate ORFs')
    update_progressbar(progressbar_file_path, 'Translate ORFs')

    if args.step_to_complete == '2':
        logger.info("Step 2 completed.")
        return

    if not args.only_calc_ogs:
        step_3_analyze_genome_completeness(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                           final_output_dir, done_files_dir, translated_orfs_dir)
        update_progressbar(progressbar_file_path, 'Calculate genomes completeness')

    if args.step_to_complete == '3':
        logger.info("Step 3 completed.")
        return

    orthogroups_file_path = step_5_full_orthogroups_inference(args, logger, times_logger, error_file_path, output_dir,
                                                              tmp_dir, done_files_dir, translated_orfs_dir,
                                                              all_proteins_fasta_path, genomes_names_path)
    update_progressbar(progressbar_file_path, 'Infer orthogroups')

    if args.step_to_complete == '5':
        logger.info("Step 5 completed.")
        return

    orphan_genes_dir = step_6_extract_orphan_genes(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                                   final_output_dir, done_files_dir, orfs_dir, orthogroups_file_path)
    update_progressbar(progressbar_file_path, 'Find orphan genes')

    if args.step_to_complete == '6':
        logger.info("Step 6 completed.")
        return

    final_orthologs_table_file_path = \
        step_7_orthologs_table_variations(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                          final_output_dir, done_files_dir, orthogroups_file_path, orphan_genes_dir)

    if args.step_to_complete == '7':
        logger.info("Step 7 completed.")
        return

    if args.only_calc_ogs:
        return

    ogs_dna_sequences_path, og_aa_sequences_path, ogs_aa_alignments_path, ogs_dna_alignments_path = \
        step_8_build_orthologous_groups_fastas(args, logger, times_logger, error_file_path,
                                               output_dir, tmp_dir, final_output_dir, done_files_dir,
                                               all_orfs_fasta_path, all_proteins_fasta_path,
                                               final_orthologs_table_file_path)
    update_progressbar(progressbar_file_path, 'Prepare orthogroups fasta files')

    if args.step_to_complete == '8':
        logger.info("Step 8 completed.")
        return

    aligned_core_proteome_reduced_file_path, core_proteome_reduced_length = step_9_extract_core_genome_and_core_proteome(
        args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir, done_files_dir,
        genomes_names_path, ogs_aa_alignments_path, ogs_dna_alignments_path)
    update_progressbar(progressbar_file_path, 'Infer core genome and proteome')

    if args.step_to_complete == '9':
        logger.info("Step 9 completed.")
        return

    if args.core_minimal_percentage == 100:
        step_10_genome_numeric_representation(args, logger, times_logger, error_file_path, output_dir, tmp_dir,
                                              final_output_dir, done_files_dir, orfs_dir, final_orthologs_table_file_path)
        update_progressbar(progressbar_file_path, 'Calculate genomes numeric representation')

    if args.step_to_complete == '10':
        logger.info("Step 10 completed.")
        return

    step_11_phylogeny(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                      done_files_dir, aligned_core_proteome_reduced_file_path, genomes_names_path, number_of_genomes,
                      core_proteome_reduced_length, data_path)
    update_progressbar(progressbar_file_path, 'Reconstruct species phylogeny')
    update_progressbar(progressbar_file_path, 'Calculate ANI (Average Nucleotide Identity)')

    if args.step_to_complete == '11':
        logger.info("Step 11 completed.")
        return

    step_12_orthogroups_annotations(args, logger, times_logger, error_file_path, output_dir, tmp_dir, final_output_dir,
                                    done_files_dir, orfs_dir, ogs_dna_sequences_path, og_aa_sequences_path,
                                    final_orthologs_table_file_path)
    update_progressbar(progressbar_file_path, 'Analyze orthogroups codon bias')
    update_progressbar(progressbar_file_path, 'Annotate orthogroups with KO terms')

    if args.step_to_complete == '12':
        logger.info("Step 12 completed.")
        return


def main(args):
    start_time = time()

    (logger, times_logger, meta_output_dir, error_file_path, progressbar_file_path, run_number, output_dir, tmp_dir,
     done_files_dir, steps_results_dir, data_path, final_output_dir_name, final_output_dir) = \
        prepare_pipeline_framework(args)

    try:
        validate_arguments(args)
        initialize_progressbar(args, progressbar_file_path)

        done_file_path = os.path.join(done_files_dir, '00_prepare_and_verify_inputs.txt')
        genomes_names_path = os.path.join(output_dir, 'genomes_names.txt')
        if not os.path.exists(done_file_path):
            prepare_and_verify_input_data(args, logger, meta_output_dir, error_file_path, data_path, genomes_names_path)
            write_to_file(logger, done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists. Skipping step...')

        update_progressbar(progressbar_file_path, 'Validate input files')

        run_main_pipeline(args, logger, times_logger, error_file_path, progressbar_file_path,
                          steps_results_dir, tmp_dir, done_files_dir,
                          data_path, genomes_names_path, final_output_dir)

        if args.step_to_complete is None and not args.only_calc_ogs and not args.do_not_copy_outputs_to_final_results_dir:
            logger.info('Zipping results folder...')
            shutil.make_archive(final_output_dir, 'zip', meta_output_dir, final_output_dir_name)
            logger.info(f'Zipped results folder to {final_output_dir}.zip')

        # remove intermediate results
        if consts.ENV == 'lsweb' and not consts.KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR:
            logger.info('Cleaning up intermediate results...')
            remove_path(logger, steps_results_dir)
            remove_path(logger, data_path)
            logger.info('Intermediate results cleaned up')

        update_progressbar(progressbar_file_path, 'Finalize results')
        state = State.Finished
    except Exception as e:
        if not os.path.exists(error_file_path):
            with open(error_file_path, 'a+') as f:
                traceback.print_exc(file=f)
        state = State.Crashed

    total_time = timedelta(seconds=int(time() - start_time))
    times_logger.info(f'Total pipeline time: {total_time}. Done')

    if consts.ENV == 'lsweb' and flask_interface_consts.SEND_EMAIL_WHEN_JOB_FINISHED_FROM_PIPELINE:
        send_email_in_pipeline_end(logger, run_number, args.email, args.job_name, state)


if __name__ == '__main__':
    print(f'sys.path is\n{sys.path}')
    arguments = get_arguments()

    main(arguments)
