import argparse
import logging
import os
import shutil
import sys
from time import time, sleep
import traceback
import json

import matplotlib.pyplot as plt
import mmap
import numpy as np
import pandas as pd

from auxiliaries.email_sender import send_email
from auxiliaries.file_writer import write_to_file
from auxiliaries.input_verifications import verify_fasta_format
from auxiliaries.pipeline_auxiliaries import measure_time, execute, wait_for_results, \
    prepare_directories, fail, submit_mini_batch, submit_batch, remove_bootstrap_values, \
    notify_admin, add_results_to_final_dir, remove_path, unpack_data, fix_illegal_chars_in_file_name, move_file
from auxiliaries.html_editor import edit_success_html, edit_failure_html, edit_progress
from auxiliaries import consts, flask_interface_consts
from auxiliaries.plots_generator import generate_violinplot, generate_bar_plot
from auxiliaries.mimic_prodigal_header import mimic_prodigal_output


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--args_json_path', help='path to a json file that contains values for arguments which will '
                                                 'override the default values. Optional.')
    parser.add_argument('--contigs_dir',
                        help='path to a folder with the genomic sequences. This folder may be zipped, as well the files'
                             ' in it.')
    parser.add_argument('--output_dir', help='relative path of directory where the output files will be written to',
                        default='outputs')
    parser.add_argument('--email', help='A notification will be sent once the pipeline is done',
                        default=consts.OWNER_EMAIL)
    parser.add_argument('--identity_cutoff', default=80,
                        help='minimum required percent of identity level (lower values will be filtered out)')
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
    parser.add_argument('--inputs_are_annotated_proteomes', action='store_true',
                        help='whether the input files are genomes or annotated proteomes')
    parser.add_argument('--qfo_benchmark', action='store_true',
                        help='whether the input files are annotated proteomes in the QfO benchmark format')
    # choices=['pupkoweb', 'pupkowebr', 'pupkolab', 'pupkolabr', 'pupkotmp', 'pupkotmpr', 'itaym', 'lilach',
    # 'bioseq', 'bental', 'oren.q', 'bioseq20.q'])
    parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                        default=consts.QUEUE_FOR_JOBS)
    parser.add_argument('--dummy_delimiter',
                        help='The queue does not "like" very long commands. A dummy delimiter is used to break each row'
                             ' into different commands of a single job',
                        default='!@#')
    parser.add_argument('--src_dir', help='source code directory',
                        default=os.path.join(consts.PROJECT_ROOT_DIR, 'pipeline'))
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

    parser.add_argument('--promoters_length', default=300,
                        help='How many basepairs upstream to the ATG will be used for the sweeps analysis',
                        type=lambda x: int(x) if int(x) >= 0 else parser.error(
                            f'Minimal number of upstream basepairs should be non-negative!'))
    parser.add_argument('--minimal_number_of_sequences_allowed_for_sweeps_analysis', default=21,
                        help='MSAs with fewer sequences will be ignore in the sweeps analysis.',
                        type=lambda x: int(x) if int(x) > 4 else parser.error(
                            f'Minimal number of sequences required for sweeps analysis is 5!'))

    args = parser.parse_args()

    # Override arguments with args_json_path content
    if args.args_json_path:
        with open(args.args_json_path, 'r') as args_json_file:
            args_json = json.load(args_json_file)
            args.__dict__.update(args_json)

    return args


def validate_arguments(args):
    if os.path.exists(args.contigs_dir):
        args.contigs_dir = args.contigs_dir.rstrip('/')
    else:
        raise ValueError(f'contigs_dir argument {args.contigs_dir} does not exist!')

    args.identity_cutoff = float(args.identity_cutoff)
    if args.identity_cutoff < 0 or args.identity_cutoff > 100:
        raise ValueError(f'identity_cutoff argument {args.identity_cutoff} has invalid value')

    args.e_value_cutoff = float(args.e_value_cutoff)
    if args.e_value_cutoff < 0 or args.e_value_cutoff > 1:
        raise ValueError(f'e_value_cutoff argument {args.e_value_cutoff} has invalid value')

    args.core_minimal_percentage = float(args.core_minimal_percentage)
    if args.core_minimal_percentage < 0 or args.core_minimal_percentage > 100:
        raise ValueError(f'core_minimal_percentage argument {args.core_minimal_percentage} has invalid value')


def prepare_pipeline_framework(args):
    meta_output_dir = os.path.join(os.path.split(args.contigs_dir)[0])

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

    run_number = os.path.join(os.path.split(meta_output_dir)[1])
    logger.info(f'run_number is {run_number}')

    if not consts.IGNORE_HTML:
        output_html_path = os.path.join(meta_output_dir, consts.RESULT_WEBPAGE_NAME)
        logger.info(f'output_html_path is {output_html_path}')

        with open(output_html_path) as f:
            html_content = f.read()
        html_content = html_content.replace('QUEUED', 'RUNNING')
        if 'progress-bar-striped active' not in html_content:
            html_content = html_content.replace('progress-bar-striped', 'progress-bar-striped active')
        with open(output_html_path, 'w') as f:
            f.write(html_content)

        output_url = os.path.join(consts.WEBSERVER_RESULTS_URL, run_number, consts.RESULT_WEBPAGE_NAME)
        logger.info(f'output_url is {output_url}')

        meta_output_url = os.path.join(consts.WEBSERVER_RESULTS_URL, run_number)
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

    return logger, times_logger, meta_output_dir, error_file_path, run_number, output_html_path, output_url, meta_output_url, \
        output_dir, tmp_dir, done_files_dir


def prepare_and_verify_input_data(args, logger, meta_output_dir, error_file_path, output_dir):
    # extract zip and detect data folder
    primary_data_path = unpack_data(logger, args.contigs_dir, meta_output_dir, error_file_path)

    for system_file in os.listdir(primary_data_path):
        if system_file.startswith(('.', '_')):
            system_file_path = os.path.join(primary_data_path, system_file)
            logger.warning(f'Removing system file: {system_file_path}')
            remove_path(logger, system_file_path)

    # copies input contigs_dir because we edit the files and want to keep the input directory as is
    data_path = os.path.join(output_dir, 'inputs')
    shutil.copytree(primary_data_path, data_path, dirs_exist_ok=True)
    logger.info(f'data_path is: {data_path}')

    # have to be AFTER system files removal (in the weird case a file name starts with a space)
    filename_prefixes = set()
    for file_name in os.listdir(data_path):

        filename_prefix = os.path.splitext(file_name)[0]
        if filename_prefix in filename_prefixes:
            error_msg = f'Two (or more) of the uploaded geonmes contain the same name (prefix), ' \
                        f'e.g., {filename_prefix}. Please make sure each file name is unique.'
            fail(logger, error_msg, error_file_path)
        filename_prefixes.add(filename_prefix)

        new_file_name = fix_illegal_chars_in_file_name(logger, file_name)
        if file_name != new_file_name:
            # illegal character in file name were found
            move_file(logger, data_path, file_name, new_file_name, error_file_path)

    number_of_genomes = len(os.listdir(data_path))
    logger.info(f'Number of genomes to analyze is {number_of_genomes}')
    logger.info(f'data_path contains the following {number_of_genomes} files: {os.listdir(data_path)}')

    # check MINimal number of genomes
    min_number_of_genomes_to_analyze = 2
    if number_of_genomes < min_number_of_genomes_to_analyze:
        error_msg = f'The dataset contains too few genomes ({consts.WEBSERVER_NAME} does comparative analysis and ' \
                    f'thus needs at least 2 genomes).'
        fail(logger, error_msg, error_file_path)

    # check MAXimal number of genomes
    max_number_of_genomes_to_analyze = 350
    if number_of_genomes > max_number_of_genomes_to_analyze and 'oren' not in args.email.lower():
        error_msg = f'The dataset contains too many genomes. {consts.WEBSERVER_NAME} allows analyzing up to ' \
                    f'{max_number_of_genomes_to_analyze} genomes due to the high resource consumption. However, ' \
                    f'upon request (and supervision), we do allow analyzing large datasets. Please contact us ' \
                    f'directly and we will do that for you.'
        fail(logger, error_msg, error_file_path)

    # must be only after the spaces removal from the species names!!
    verification_error = verify_fasta_format(logger, data_path)
    if verification_error:
        remove_path(logger, data_path)
        fail(logger, verification_error, error_file_path)

    return data_path, number_of_genomes


def run_main_pipeline(args, logger, times_logger, meta_output_dir, error_file_path, run_number, output_html_path,
                      output_dir, tmp_dir, done_files_dir, data_path, number_of_genomes):
    # 1.	drop_plasmids.py
    # Input: (1) an input path for a fasta file with contigs/full genome
    # Output: overrides the original file without plasmids
    # Can be parallelized on cluster
    if args.filter_out_plasmids:
        step_number = '01'
        logger.info(f'Step {step_number}: {"_" * 100}')
        step_name = f'{step_number}_drop_plasmids'
        script_path = os.path.join(args.src_dir, 'steps/drop_plasmids.py')
        filtered_inputs_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
        done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
        if not os.path.exists(done_file_path):
            logger.info('Filtering plasmids out...')
            all_cmds_params = []
            for fasta_file in os.listdir(data_path):
                single_cmd_params = [os.path.join(data_path, fasta_file),
                                     os.path.join(filtered_inputs_dir, fasta_file)]
                all_cmds_params.append(single_cmd_params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                       num_of_cmds_per_job=5,
                                                       job_name_suffix='drop_plasmids',
                                                       queue_name=args.queue_name)

            wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                             num_of_batches, error_file_path, email=args.email)
            write_to_file(logger, done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists. Skipping step...')
        data_path = filtered_inputs_dir
        edit_progress(output_html_path, progress=5)

    # 2.	search_orfs.py
    # Input: (1) an input path for a fasta file with contigs/full genome (2) an output file path
    #                (with a suffix as follows: i_genes.fasta. especially relevant for the wrapper).
    # Output: a fasta file where under each header, there’s a single gene.
    # Can be parallelized on cluster
    # Prodigal ignores newlines in the middle of a sequence, namely, >bac1\nAAA\nAA\nTTT >bac1\nAAAAATTT will be
    # analyzed identically.
    step_number = orfs_step_number = '02'  # using orfs_step_number for downstream analysis (extract promoter and gene sequences)
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_ORFs'
    script_path = os.path.join(args.src_dir, 'steps/search_orfs.py')
    orfs_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        if not args.inputs_are_annotated_proteomes:
            logger.info('Extracting ORFs...')
            all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
            for fasta_file in os.listdir(data_path):
                fasta_file_prefix = os.path.splitext(fasta_file)[0]
                output_file_name = f'{fasta_file_prefix}.{step_name}'

                single_cmd_params = [f'"{os.path.join(data_path, fasta_file)}"',
                                     os.path.join(orfs_dir, output_file_name)]
                all_cmds_params.append(single_cmd_params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                       num_of_cmds_per_job=5,
                                                       job_name_suffix='search_orfs',
                                                       queue_name=args.queue_name,
                                                       required_modules_as_list=[consts.PRODIGAL])

            wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                             num_of_batches, error_file_path, email=args.email)
        else:  # inputs are annotated proteomes
            logger.info('Inputs are already annotated proteomes. Skipping step 02.')
            shutil.copytree(data_path, orfs_dir, dirs_exist_ok=True)
            mimic_prodigal_output(orfs_dir, step_name)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=10)

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
                error_msg = f'{consts.WEBSERVER_NAME} could not detect any ORFs in:\n<br> '
                missing_orfs = 1
            # add genome-without-orfs name
            error_msg += f'{os.path.splitext(file)[0]}\n<br>'

    if missing_orfs:
        error_msg += f'\n<br>Please remove the abovementioned from the dataset and re-submit your job. In general, ' \
                     f'it is recommended to use genomic files that contain at least 20K base pairs (each).'
        fail(logger, error_msg, error_file_path)

    # 3.  create_mmseqs2_DB.py
    # Input: path to gene file to create DB from
    # Output: DB of the input file for mmseqs2
    # Can be parallelized on cluster
    step_number = '03'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = orfs_dir
    step_name = f'{step_number}_dbs'
    script_path = os.path.join(args.src_dir, 'steps/create_mmseqs2_DB.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Creating DBs...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for fasta_file in os.listdir(previous_pipeline_step_output_dir):
            file_path = os.path.join(previous_pipeline_step_output_dir, fasta_file)
            fasta_file_prefix = os.path.splitext(fasta_file)[0]
            output_file_name = f'{fasta_file_prefix}'
            output_prefix = os.path.join(pipeline_step_output_dir, output_file_name)

            single_cmd_params = [file_path,
                                 output_prefix,
                                 output_prefix,  # instead of tmp_dir
                                 '-t']  # translate to peptides # Should we let the user decide?
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=10,
                                                   job_name_suffix='DB_creation',
                                                   queue_name=consts.QUEUE_FOR_MMSEQS_COMMANDS,
                                                   required_modules_as_list=[consts.MMSEQS])

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=15)

    # send a subjob that removes all mmseqs intermediate *FOLDERS* (e.g., 3465136234521948 etc..) in tmp_dir
    # does not remove (sge/cmds/log) files. only folders.
    # MMSEQS generates tons of intermediate files that abuse the inodes so they are deleted once the step is over!
    submit_mini_batch(logger, os.path.join(args.src_dir, 'auxiliaries/remove_tmp_folders.py'),
                      [[pipeline_step_tmp_dir]], pipeline_step_tmp_dir, args.queue_name,
                      job_name='remove_dirs_from_tmp')

    # 4.	mmseqs2_all_vs_all.py
    # Input: (1) 2 input paths for 2 (different) genome files (query and target), g1 and g2
    #        (2) an output file path (with a suffix as follows: i_vs_j.tsv. especially relevant for the wrapper).
    # Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
    # Precisely, for each gene x in g1, query x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
    # Can be parallelized on cluster
    step_number = '04'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_all_vs_all_analysis'
    script_path = os.path.join(args.src_dir, 'steps/mmseqs2_all_vs_all.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info(f'Querying all VS all (using mmseqs2)...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for db1 in os.listdir(orfs_dir):
            strain1_name = os.path.splitext(db1)[0]
            for db2 in os.listdir(orfs_dir):
                strain2_name = os.path.splitext(db2)[0]

                logger.debug(f'{"#" * 100}\nGot {strain1_name} and {strain2_name}')
                if strain1_name >= strain2_name:
                    logger.debug(f'Skipping {strain1_name} >= {strain2_name}')
                    continue  # no need to query strain against itself or a pair that was already seen
                logger.info(f'Querying {strain1_name} against {strain2_name}')

                query_aa_db = os.path.splitext(os.path.join(previous_pipeline_step_output_dir, db1))[0] + '_aa'

                target_aa_db = os.path.splitext(os.path.join(previous_pipeline_step_output_dir, db2))[0] + '_aa'

                aln_offsetted_db = os.path.join(pipeline_step_tmp_dir,
                                                f'{strain1_name}_vs_{strain2_name}.alnOffsettedDB')

                output_file_name = f'{strain1_name}_vs_{strain2_name}.m8'
                output_file_path = os.path.join(pipeline_step_output_dir, output_file_name)

                single_cmd_params = [query_aa_db,
                                     target_aa_db,
                                     aln_offsetted_db,
                                     f'{pipeline_step_tmp_dir}/{strain1_name}_vs_{strain2_name}',
                                     output_file_path]

                all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100 if len(os.listdir(orfs_dir)) > 25 else 5,
                                                   job_name_suffix='rbh_analysis',
                                                   queue_name=consts.QUEUE_FOR_MMSEQS_COMMANDS,
                                                   required_modules_as_list=[consts.MMSEQS])

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=20)

    # send a subjob that removes (now instead of at the end) all mmseqs intermediate *FOLDERS*
    # (e.g., 3465136234521948 etc..) in tmp_dir
    # does not remove (sge/cmds/log) files. only folders.
    submit_mini_batch(logger, os.path.join(args.src_dir, 'auxiliaries/remove_tmp_folders.py'),
                      [[previous_pipeline_step_output_dir]],
                      pipeline_step_tmp_dir, args.queue_name, job_name='remove_m8_files')

    # 5.	filter_rbh_results.py
    # Input: (1) a path for a i_vs_j.tsv file (2) an output path (with a suffix as follows: i_vs_j_filtered.tsv.
    #            especially relevant for the wrapper).
    # Output: the same format of the input file containing only pairs that passed the filtration. For each row in
    #         the input file (pair of genes), apply the following filters:
    # 1. at least X% similarity
    # 2. at least X% of the length
    # 3.# write each pair to the output file if it passed all the above filters.
    # Can be parallelized on cluster
    step_number = '05'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_blast_filtered'
    script_path = os.path.join(args.src_dir, 'steps/filter_rbh_results.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Filtering all vs all results...\n')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for blast_results_file in os.listdir(previous_pipeline_step_output_dir):
            fasta_file_prefix = os.path.splitext(blast_results_file)[0]
            output_file_name = f'{fasta_file_prefix}.{step_name}'

            single_cmd_params = [os.path.join(previous_pipeline_step_output_dir, blast_results_file),
                                 os.path.join(pipeline_step_output_dir, output_file_name),
                                 f'--identity_cutoff {args.identity_cutoff / 100}',
                                 # needs to be normalized between 0 and 1
                                 f'--e_value_cutoff {args.e_value_cutoff}']
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100 if len(os.listdir(orfs_dir)) > 100 else 50,
                                                   job_name_suffix='hits_filtration',
                                                   queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=25)

    # 6. concatenate_reciprocal_hits
    # Input: path to folder with all reciprocal hits files
    # Output: concatenated file of all reciprocal hits files
    # CANNOT be parallelized on cluster
    step_number = '06'
    logger.info(f'Step {step_number}: {"_" * 100}')
    all_reciprocal_hits_file = os.path.join(output_dir, 'concatenated_all_reciprocal_hits.txt')
    if not os.path.exists(all_reciprocal_hits_file):
        logger.info('Concatenating reciprocal hits...')
        for filtered_hits_file in os.listdir(pipeline_step_output_dir):
            execute(logger, f'cat {pipeline_step_output_dir}/{filtered_hits_file} >> {all_reciprocal_hits_file}',
                    process_is_string=True)
        # avoid cat {pipeline_step_output_dir}/* because arguments list might be too long!
        # No need to wait...
    else:
        logger.info(f'done file {all_reciprocal_hits_file} already exists. Skipping step...')
    edit_progress(output_html_path, progress=30)

    # 6b   extract_orphan_genes.py
    step_number = '6b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orphan_genes'
    script_path = os.path.join(args.src_dir, 'steps/extract_orphan_genes.py')
    orphan_genes_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting orphan genes...')
        job_name = os.path.split(script_path)[-1]
        params = [orfs_dir,
                  all_reciprocal_hits_file,
                  orphan_genes_dir]
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, args.queue_name, job_name=job_name)
        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=35)

    # 7.	construct_putative_orthologs_table.py
    # Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
    # Output: updates the table with the info from the reciprocal hit file.
    # CANNOT be parallelized on cluster
    step_number = '07'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_putative_table'
    script_path = os.path.join(args.src_dir, 'steps/construct_putative_orthologs_table.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    putative_orthologs_table_path = os.path.join(pipeline_step_output_dir, 'putative_orthologs_table.txt')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing putative orthologs table...')
        job_name = os.path.split(script_path)[-1]
        params = [all_reciprocal_hits_file,
                  putative_orthologs_table_path]
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, args.queue_name, job_name=job_name)
        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=35)

    # 8   prepare_files_for_mcl.py
    # Input: (1) a path for a concatenated all reciprocal hits file (2) a path for a putative orthologs file (3) a path for an output folder
    # Output: an input file for MCL for each putative orthologs group
    # CANNOT be parallelized on cluster (if running on the concatenated file)
    step_number = '08'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_mcl_input_files'
    script_path = os.path.join(args.src_dir, 'steps/prepare_files_for_mcl.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Preparing files for MCL...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        clusters_to_prepare_per_job = 10
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

            single_cmd_params = [all_reciprocal_hits_file,
                                 putative_orthologs_table_path,
                                 first_mcl,
                                 last_mcl,
                                 pipeline_step_output_dir]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=5,
                                                   # *times* the number of clusters_to_prepare_per_job above. 50 in total per batch!
                                                   job_name_suffix='mcl_preparation',
                                                   queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=40)

    # 9.	run_mcl.py
    # Input: (1) a path to an MCL input file (2) a path to MCL's output.
    # Output: MCL analysis.
    # Can be parallelized on cluster
    step_number = '09'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_mcl_analysis'
    script_path = os.path.join(args.src_dir, 'steps/run_mcl.py')
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

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100,
                                                   job_name_suffix='mcl_execution',
                                                   queue_name=args.queue_name,
                                                   required_modules_as_list=[consts.MCL])

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=45)

    # 10.	verify_cluster.py
    # Input: (1) mcl analysis file (2) a path to which the file will be moved if relevant (3) optional: maximum number of clusters allowed [default=1]
    # Output: filter irrelevant clusters by moving the relevant to an output directory
    # Can be parallelized on cluster
    step_number = '10'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_verified_clusters'
    script_path = os.path.join(args.src_dir, 'steps/verify_cluster.py')
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Verifying clusters...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for putative_orthologs_group in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_group_prefix = os.path.splitext(putative_orthologs_group)[0]
            job_name = os.path.split(putative_orthologs_group_prefix)[-1]

            single_cmd_params = [f'"{os.path.join(previous_pipeline_step_output_dir, putative_orthologs_group)}"',
                                 f'"{os.path.join(pipeline_step_output_dir, putative_orthologs_group)}"']
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100,
                                                   job_name_suffix='clusters_verification',
                                                   queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=50)

    logger.info(f'A total of {len(os.listdir(previous_pipeline_step_output_dir))} clusters were verified.')
    logger.debug(f'The verified clusters in {previous_pipeline_step_output_dir} are the following:')
    logger.debug(os.listdir(previous_pipeline_step_output_dir))

    # 11.	construct_final_orthologs_table.py
    # Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
    # Output: aggregates all the well-clustered OGs to the final table.
    step_number = '11'
    logger.info(f'Step {step_number}: {"_" * 100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    step_name = f'{step_number}_final_table'
    script_path = os.path.join(args.src_dir, 'steps/construct_final_orthologs_table.py')
    final_orthologs_table_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    final_orthologs_table_file_path = os.path.join(final_orthologs_table_path, 'final_orthologs_table.csv')
    phyletic_patterns_path = os.path.join(final_orthologs_table_path, 'phyletic_pattern.fas')
    orthoxml_path = os.path.join(final_orthologs_table_path, 'orthologs.orthoxml')
    ortholog_pairs_path = os.path.join(final_orthologs_table_path, 'ortholog_pairs.tsv')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing final orthologs table...')
        params = [putative_orthologs_table_path,
                  previous_pipeline_step_output_dir,
                  final_orthologs_table_file_path,
                  phyletic_patterns_path,
                  orthoxml_path,
                  ortholog_pairs_path]
        if args.qfo_benchmark:
            params += ['--qfo_benchmark']
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, args.queue_name, job_name='final_ortholog_groups')
        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path,
                         error_message='No ortholog groups were detected in your dataset. Please try to lower '
                                       'the similarity parameters (see Advanced Options in the submission page) '
                                       'and re-submit your job.', email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=55)

    # extract orthologs table header for sequence extraction later on
    with open(final_orthologs_table_file_path) as f:
        header_line = f.readline()
        first_delimiter_index = header_line.index(consts.CSV_DELIMITER)
        final_table_header = header_line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"

    # extract of strains names for core genome analysis later on
    strains_names_path = os.path.join(final_orthologs_table_path, 'strains_names.txt')
    strains_names = final_table_header.split(consts.CSV_DELIMITER)
    with open(strains_names_path, 'w') as f:
        f.write('\n'.join(strains_names) + '\n')

    # extract number of strains for core genome analysis later on
    num_of_strains_path = os.path.join(final_orthologs_table_path, 'num_of_strains.txt')
    num_of_strains = len(strains_names)
    with open(num_of_strains_path, 'w') as f:
        f.write(f'{num_of_strains}\n')

    # 12.	extract_orfs.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    step_number = '12'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthologs_groups_dna_sequences'
    script_path = os.path.join(args.src_dir, 'steps/extract_orfs.py')
    orthologs_dna_sequences_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        # og_number = 0
        with open(final_orthologs_table_file_path) as f:
            f.readline()  # skip header
            for line in f:
                first_delimiter_index = line.index(consts.CSV_DELIMITER)
                og_name = line[:first_delimiter_index]
                cluster_members = line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"
                output_file_name = og_name

                single_cmd_params = [orfs_dir,
                                     f'"{final_table_header}"',
                                     # should be flanked by quotes because it might contain spaces...
                                     f'"{cluster_members}"',
                                     # should be flanked by quotes because it might contain spaces...
                                     f'"{og_name}"',  # should be flanked by quotes because it might contain spaces...
                                     os.path.join(orthologs_dna_sequences_dir_path, f'{output_file_name}_dna.fas')]
                all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=250,
                                                   job_name_suffix='orfs_extraction',
                                                   queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=60)

    # 12_5  codon_bias.py
    # Input: ORF dir and OG dir
    # Output: W_vectors, CAI for each OG
    step_number = '12_5'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_codon_bias'
    script_path = os.path.join(args.src_dir, 'tests/codon_bias/main.py')
    codon_bias_output_dir_path, codon_bias_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Analyzing codon bias...')
        params = [
            orfs_dir,
            orthologs_dna_sequences_dir_path,
            codon_bias_output_dir_path,
            codon_bias_tmp_dir,
            args.src_dir,
            args.queue_name,
            error_file_path
        ]
        submit_mini_batch(logger, script_path, [params], codon_bias_tmp_dir,
                          args.queue_name, job_name='codon_bias')

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], codon_bias_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=65)

    # 13.  translate_fna_to_faa.py
    # Input: path to fna file and an faa file
    # Output: translate the fna to protein and write to the faa file
    # Can be parallelized on cluster
    step_number = '13'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthologs_groups_aa_sequences'
    script_path = os.path.join(args.src_dir, 'steps/translate_fna_to_faa.py')
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

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=250,
                                                   job_name_suffix='dna_translation',
                                                   queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=65)

    # 14.	align_orthologs_group.py
    # Input: (1) A path to an unaligned amino acid sequences file (2) An output file path
    # Output: aligned sequences
    # Can be parallelized on cluster
    step_number = '14'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_aligned_aa_orthologs_groups'
    script_path = os.path.join(args.src_dir, 'steps/align_orthologs_group.py')
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
                                                   pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100,
                                                   job_name_suffix='genes_alignment',
                                                   queue_name=args.queue_name,
                                                   required_modules_as_list=[consts.MAFFT])

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=70)

    # 15.	extract aligned_core_proteome.py
    step_number = '15'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_aligned_core_proteome'
    script_path = os.path.join(args.src_dir, 'steps/extract_core_genome.py')
    aligned_core_proteome_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    aligned_core_proteome_file_path = os.path.join(aligned_core_proteome_path, 'aligned_core_proteome.fas')
    core_ogs_names_file_path = os.path.join(aligned_core_proteome_path, 'core_ortholog_groups_names.txt')
    core_length_file_path = os.path.join(aligned_core_proteome_path, 'core_length.txt')
    number_of_core_members = os.path.join(aligned_core_proteome_path, 'number_of_core_genes.txt')
    with open(num_of_strains_path) as f:
        num_of_strains = f.read().rstrip()
    if not os.path.exists(done_file_path):
        logger.info('Extracting aligned core proteome...')

        params = [aa_alignments_path, num_of_strains, strains_names_path,
                  aligned_core_proteome_file_path,
                  core_ogs_names_file_path,
                  core_length_file_path,
                  number_of_core_members,
                  f'--core_minimal_percentage {args.core_minimal_percentage}']  # how many members induce a core group?
        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir,
                          args.queue_name, job_name='core_proteome')

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path, email=args.email)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=75)

    # 16.	reconstruct_species_phylogeny.py
    step_number = '16'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_species_phylogeny'
    script_path = os.path.join(args.src_dir, 'steps/reconstruct_species_phylogeny.py')
    phylogeny_path, phylogeny_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    phylogenetic_raw_tree_path = os.path.join(phylogeny_path, 'final_species_tree.txt')
    start_tree = time()
    if not os.path.exists(phylogenetic_raw_tree_path):
        logger.info('Reconstructing species phylogeny...')

        if 'oren' in args.email or 'pupkoweb' not in args.queue_name:
            num_of_cpus = 20  # anyhow it can't be more than 20! o.w., "qsub: Job violates queue and/or server resource limits"
        else:
            # better not to use more than 10 cpus, routinely, when running on pupkoweb
            num_of_cpus = 10
        params = [aligned_core_proteome_file_path,
                  phylogenetic_raw_tree_path,
                  f'--cpu {num_of_cpus}']
        if args.outgroup:
            if args.outgroup in strains_names:
                params += [f'--outgroup {args.outgroup}']
            else:
                logger.info(f'Outgroup {args.outgroup} was specified but it is not one of the input species:\n'
                            f'{",".join(sorted(strains_names))}\nAn unrooted tree is going to be reconstructed')
        if args.bootstrap == 'yes' and (
                number_of_genomes < 150 or 'oren' in args.email):  # allow bootstrap only for less than 150 genomes
            params += ['--num_of_bootstrap_iterations 100']

        submit_mini_batch(logger, script_path, [params], phylogeny_tmp_dir,
                          args.queue_name, job_name='tree_reconstruction',
                          required_modules_as_list=[consts.RAXML], num_of_cpus=num_of_cpus)
        # no need to wait now. Wait before plotting the tree!
        # Still, in order to get more than one thread, allow few seconds to the subjobs to be submitted
        # (before the next batch is submitted and pounds the cluster)
        sleep(60)
    else:
        logger.info(f'Raw tree file {phylogenetic_raw_tree_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=80)

    # 16b.	genome_numeric_representation.py
    step_number = '16b'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_genome_numeric_representation'
    script_path = os.path.join(args.src_dir, 'steps/genome_numeric_representation.py')
    numeric_representation_output_dir, numeric_representation_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    core_genome_numeric_representation_file_path = os.path.join(numeric_representation_output_dir, 'core_genome_numeric_representation.txt')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        params = [final_orthologs_table_file_path,
                  orfs_dir,
                  core_genome_numeric_representation_file_path]
        submit_mini_batch(logger, script_path, [params], numeric_representation_tmp_dir,
                          args.queue_name, job_name='numeric_representation')

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], numeric_representation_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path, email=args.email)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=90)

    # 17.	extract_orfs_statistics.py
    # Input: (1) A path to ORFs file (2) An output path to ORFs counts (3) An output path to GC content
    # Output: write the number of ORFs and GC content to the output files (respectively)
    # Can be parallelized on cluster
    step_number = '17'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orfs_statistics'
    script_path = os.path.join(args.src_dir, 'steps/extract_orfs_statistics.py')
    orfs_statistics_path, orfs_statistics_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    start_orf_stats = time()
    if not os.path.exists(done_file_path):
        logger.info('Collecting orfs counts...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for file in os.listdir(orfs_dir):
            orf_path = os.path.join(orfs_dir, file)
            strain_name = os.path.splitext(file)[0]
            orfs_count_output_path = os.path.join(orfs_statistics_path, f'{strain_name}.orfs_count')
            gc_content_output_path = os.path.join(orfs_statistics_path, f'{strain_name}.gc_content')

            single_cmd_params = [orf_path, orfs_count_output_path, gc_content_output_path]
            all_cmds_params.append(single_cmd_params)

        num_of_expected_orfs_results, example_cmd = submit_batch(logger, script_path, all_cmds_params,
                                                                 orfs_statistics_tmp_dir,
                                                                 num_of_cmds_per_job=200,
                                                                 job_name_suffix='orfs_statistics',
                                                                 queue_name=args.queue_name)

        # no need to wait now. Wait before plotting the statistics!
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 18.	induce_dna_msa_by_aa_msa.py
    # Input: (1) An aa alignment (2) An unaligned dna file (3) An output file path
    # Output: write the codon alignment induced by the aa alignment to the output file
    # Can be parallelized on cluster
    step_number = '18'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_induce_dna_msa_by_aa_msa'
    script_path = os.path.join(args.src_dir, 'steps/induce_dna_msa_by_aa_msa.py')
    num_of_expected_induced_results = 0
    dna_alignments_path, induced_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    start_induced = time()
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
                                                                    induced_tmp_dir,
                                                                    num_of_cmds_per_job=250,
                                                                    job_name_suffix='induce_msa',
                                                                    queue_name=args.queue_name)

        # no need to wait now. Wait before moving the results dir!
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 19.	extract_groups_sizes_frequency
    step_number = '19'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_groups_sizes_frequency'
    group_sizes_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    groups_sizes_frequency_file_prefix = os.path.join(group_sizes_path, 'groups_sizes_frequency')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Collecting sizes...')
        group_sizes = []
        with open(final_orthologs_table_file_path) as f:
            final_table_header = f.readline().rstrip()
            for line in f:
                cluster_members = line.rstrip()
                size = sum(
                    bool(item) for item in cluster_members.split(consts.CSV_DELIMITER))  # count non empty entries
                size -= 1  # don't count group name
                group_sizes.append(str(size))

        groups_sizes_frequency_raw_file_path = groups_sizes_frequency_file_prefix + '.txt'
        with open(groups_sizes_frequency_raw_file_path, 'w') as f:
            f.write('\n'.join(group_sizes))

        groups_sizes_frequency_png_file_path = groups_sizes_frequency_file_prefix + '.png'
        generate_bar_plot(groups_sizes_frequency_raw_file_path, groups_sizes_frequency_png_file_path,
                          xlabel='Orthologous group size', ylabel='Count')

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 20.	plot_orfs_statistics
    step_number = '20'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orfs_plots'
    orfs_plots_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orfs_counts_frequency_file = os.path.join(orfs_plots_path, 'orfs_counts.txt')
    orfs_counts_per_genome_file = os.path.join(orfs_plots_path, 'orfs_counts.json')
    orfs_gc_content_file = os.path.join(orfs_plots_path, 'orfs_gc_contents.txt')
    orfs_gc_content_per_genome_file = os.path.join(orfs_plots_path, 'orfs_gc_contents.json')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        wait_for_results(logger, times_logger, 'extract_orfs_statistics.py', orfs_statistics_tmp_dir,
                         num_of_expected_orfs_results, error_file_path=error_file_path, start=start_orf_stats,
                         email=args.email)

        logger.info('Concatenating orfs counts and gc contents...')

        gc_content = {}
        orf_count = {}
        for file_name in os.listdir(orfs_statistics_path):
            strain_name = file_name.split('.')[0]
            file_type = file_name.split('.')[1]
            file_path = os.path.join(orfs_statistics_path, file_name)
            with open(file_path, 'r') as f:
                content = f.read()
            if file_type == "gc_content":
                gc_content[strain_name] = float(content)
            if file_type == "orfs_count":
                orf_count[strain_name] = int(content)

        with open(orfs_gc_content_per_genome_file, 'w') as f:
            json.dump(gc_content, f)
        with open(orfs_counts_per_genome_file, 'w') as f:
            json.dump(orf_count, f)

        execute(logger, f'cat {orfs_statistics_path}/*.orfs_count > {orfs_counts_frequency_file}', process_is_string=True)
        execute(logger, f'cat {orfs_statistics_path}/*.gc_content > {orfs_gc_content_file}', process_is_string=True)

        # No need to wait...
        logger.info('Plotting violins...')
        orfs_counts_frequency_png_file_path = orfs_counts_frequency_file.replace('txt', 'png')
        generate_violinplot(orfs_counts_frequency_file, orfs_counts_frequency_png_file_path,
                         xlabel='\nORF count per genome', dpi=300)

        orfs_gc_content_png_file_path = orfs_gc_content_file.replace('txt', 'png')
        generate_violinplot(orfs_gc_content_file, orfs_gc_content_png_file_path,
                         xlabel='\nGC content per genome', dpi=300)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')
    edit_progress(output_html_path, progress=85)

    # Final step_number: gather relevant results, zip them together and update html file
    logger.info(f'FINAL STEP: {"_" * 100}')
    final_output_dir_name = f'{consts.WEBSERVER_NAME}_outputs'
    final_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, final_output_dir_name)
    done_file_path = os.path.join(done_files_dir, f'{final_output_dir_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Gathering results to final output dir...')

        # move ORFs folder
        add_results_to_final_dir(logger, orfs_dir, final_output_dir, copy=True)

        # move orphan genes
        add_results_to_final_dir(logger, orphan_genes_dir, final_output_dir, copy=True)

        # move orthologs table
        add_results_to_final_dir(logger, final_orthologs_table_path, final_output_dir, copy=True)

        # move unaligned dna sequences
        add_results_to_final_dir(logger, orthologs_dna_sequences_dir_path, final_output_dir)

        # move unaligned aa sequences
        add_results_to_final_dir(logger, orthologs_aa_sequences_dir_path, final_output_dir)

        # move aligned aa sequences
        add_results_to_final_dir(logger, aa_alignments_path, final_output_dir)

        # move core proteome dir
        add_results_to_final_dir(logger, aligned_core_proteome_path, final_output_dir)

        # move groups sizes
        add_results_to_final_dir(logger, group_sizes_path, final_output_dir)

        # move orfs statistics
        add_results_to_final_dir(logger, orfs_statistics_path, final_output_dir)

        # move orfs plot
        add_results_to_final_dir(logger, orfs_plots_path, final_output_dir)

        if num_of_expected_induced_results > 0:  # can be 0 when re-running
            wait_for_results(logger, times_logger, 'induce_dna_msa_by_aa_msa.py', induced_tmp_dir,
                             num_of_expected_results=num_of_expected_induced_results,
                             error_file_path=error_file_path, start=start_induced, email=args.email)
        # move induced dna sequences
        add_results_to_final_dir(logger, dna_alignments_path, final_output_dir)

        # wait for the raw tree here
        wait_for_results(logger, times_logger, 'reconstruct_species_phylogeny.py', phylogeny_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path,
                         start=start_tree, email=args.email)
        # move species tree dir
        add_results_to_final_dir(logger, phylogeny_path, final_output_dir)

        # move numeric representation dir
        add_results_to_final_dir(logger, numeric_representation_output_dir, final_output_dir)

        edit_progress(output_html_path, progress=98)

        logger.info('Zipping results folder...')
        shutil.make_archive(final_output_dir, 'zip', final_output_dir)

        logger.info(f'Moving results to parent dir... ({meta_output_dir})')
        try:
            shutil.move(f'{final_output_dir}.zip', meta_output_dir)
        except shutil.Error as e:
            logger.error(e.args[0])
        try:
            shutil.move(final_output_dir, meta_output_dir)
        except shutil.Error as e:
            logger.error(e.args[0])

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    logger.info('Editing results html...')
    edit_success_html(logger, output_html_path, meta_output_dir, final_output_dir_name, run_number)

    edit_progress(output_html_path, progress=100, active=False)

    return orfs_dir, orfs_step_number, final_orthologs_table_file_path, phylogenetic_raw_tree_path, \
        final_output_dir_name


def report_error_in_main_pipeline_to_admin(logger, e, meta_output_dir, error_file_path, run_number, output_html_path, meta_output_url):
    error_msg = f'{consts.WEBSERVER_NAME} failed :('
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

    notify_admin(meta_output_dir, meta_output_url, run_number)


def report_main_pipeline_result_to_user(args, logger, status, total_time, output_url, run_number):
    msg = f'M1CR0B1AL1Z3R pipeline {status}'
    if status == 'is done':
        msg += f' (Took {measure_time(total_time)}).\nResults can be found at {output_url}.\nPlease note that the ' \
               f'results will be kept in the server for three months.'
    else:
        msg += f'. For further information please visit: {output_url}'
    logger.info(msg)

    logger.info(f'Sending a notification email to {args.email}')
    try:
        send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>', args.email,
                   subject=f'{consts.WEBSERVER_NAME} run number {run_number} {status}.', content=msg)
    except:
        logger.error(f'\nFailed sending notification to {args.email}\n')


def run_pipeline_extensions(args, logger, times_logger, error_file_path, run_number, output_dir, tmp_dir, done_files_dir, data_path, orfs_dir,
                            orfs_step_number, final_orthologs_table_file_path, phylogenetic_raw_tree_path,
                            final_output_dir_name):
    logger.info('\n\n\n')
    logger.info('#' * 100)
    logger.info('Executing additional modules for Sweeps analysis...')
    # 21.	extract_promoters_and_orfs
    # Input: (1) A path to a genome (2) Prodigal output file with ORFs coordinates
    # Output: A fasta file containing promoters+genes (all ORFS of the given coordinates + k[=300] upstream bases)
    # Can be parallelized on cluster
    step_number = '21'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_extract_promoters_and_orfs'
    script_path = os.path.join(args.src_dir, 'steps/extract_promoters_and_orfs.py')
    pipeline_step_output_dir_21, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting promoters...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for genome_file in os.listdir(data_path):
            genome_path = os.path.join(data_path, genome_file)
            genome_file_prefix = os.path.splitext(genome_file)[0]
            orfs_path = os.path.join(orfs_dir, f'{genome_file_prefix}.{orfs_step_number}_ORFs')
            output_prefix = os.path.join(pipeline_step_output_dir_21, f'{genome_file_prefix}.promoter_and_orf')

            single_cmd_params = [genome_path, orfs_path, output_prefix,
                                 f'--promoters_length {args.promoters_length}']  # (optional)
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=5,
                                                   job_name_suffix='promoter_gene_alignment',
                                                   queue_name=args.queue_name,
                                                   required_modules_as_list=[consts.MAFFT])

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 22.	extract_orfs.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    # Can be parallelized on cluster
    step_number = '22'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthologs_groups_dna_sequences'
    script_path = os.path.join(args.src_dir, 'steps/extract_orfs.py')
    pipeline_step_output_dir_22, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    assert os.path.exists(
        pipeline_step_output_dir_22), f'Failed to create output folder. {pipeline_step_output_dir_22} does not exist!'
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        with open(final_orthologs_table_file_path) as f:
            logger.info('Extracting orthologs groups sequences according to final orthologs table...')
            all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
            header_line = f.readline()
            first_delimiter_index = header_line.index(consts.CSV_DELIMITER)
            final_table_header = header_line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"
            for line in f:
                first_delimiter_index = line.index(consts.CSV_DELIMITER)
                og_name = line[:first_delimiter_index]
                cluster_members = line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"

                number_of_sequences_in_cluster = sum(
                    1 for member in cluster_members.split(consts.CSV_DELIMITER) if member)
                if number_of_sequences_in_cluster < args.minimal_number_of_sequences_allowed_for_sweeps_analysis:
                    # ignore msa that contains less than the minimal allowed sequences for sweeps analysis
                    logger.info(
                        f'Orthologs group #{og_name} contains {number_of_sequences_in_cluster} members and thus '
                        f'will not be included in the sweeps analysis (at least '
                        f'{args.minimal_number_of_sequences_allowed_for_sweeps_analysis} '
                        f'members are required)')
                    continue

                logger.debug(f'Orthologs group #{og_name} contains {number_of_sequences_in_cluster} sequences and'
                             f' thus will be included in the sweeps analysis')
                og_sequences_path = os.path.join(pipeline_step_output_dir_22, f'{og_name}_dna.fas')

                single_cmd_params = [pipeline_step_output_dir_21,
                                     f'"{final_table_header}"',
                                     # should be flanked by quotes because it might contain spaces...
                                     f'"{cluster_members}"',
                                     # should be flanked by quotes because it might contain spaces...
                                     f'"{og_name}"',
                                     # should be flanked by quotes because it might contain spaces...
                                     og_sequences_path]
                all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100,
                                                   job_name_suffix='extract_sequences',
                                                   queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 23.	align_orthologs_group.py
    # Input: (1) A path to an unaligned sequences file (2) An output file path
    # Output: aligned sequences
    step_number = '23'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_aligned_dna_orthologs_groups_with_promoter'
    script_path = os.path.join(args.src_dir, 'steps/align_orthologs_group.py')
    pipeline_step_output_dir_23, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Aligning orthologs groups...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for og_file in os.listdir(pipeline_step_output_dir_22):
            og_path = os.path.join(pipeline_step_output_dir_22, og_file)
            alignment_path = os.path.join(pipeline_step_output_dir_23, og_file.replace('.fas', '_mafft.fas'))

            single_cmd_params = [og_path, alignment_path, '--type nuc']
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100,
                                                   job_name_suffix='align_promoter_and_gene',
                                                   queue_name=args.queue_name,
                                                   required_modules_as_list=[consts.MAFFT])

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 24.	adjust_tree_to_msa.py
    # Input: (1) A path to an MSA file to which the tree should be adjusted
    #        (2) A path to a species tree that contains (at least) all the species in the input MSA
    #        (3) A path to a folder in which a txt file (with the same name as the msa_file) will be created. All the species that do not appear in the msa (and thus will be removed) will be written to the file that was created',
    #        (4) A path to a file in which the pruned tree will be written
    # Output: an adjusted tree per msa at (4)
    step_number = '24'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_prunned_trees'
    script_path = f'/bioseq/sincopa/adjust_tree_to_msa.py'
    pipeline_step_output_dir_24, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Removing bootstrap values from tree...')
        phylogenetic_raw_tree_path = phylogenetic_raw_tree_path.replace('outputs',
                                                                        final_output_dir_name)  # the tree was moved to the final dir..
        phylogenetic_raw_tree_path_without_bootstrap_values = phylogenetic_raw_tree_path.replace('.txt',
                                                                                                 '_no_bootstrap.txt')
        remove_bootstrap_values(phylogenetic_raw_tree_path, phylogenetic_raw_tree_path_without_bootstrap_values)

        logger.info('Pruning trees...')
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for msa_file in os.listdir(pipeline_step_output_dir_23):
            msa_path = os.path.join(pipeline_step_output_dir_23, msa_file)
            pruned_tree_path = os.path.join(pipeline_step_output_dir_24, msa_file.replace('.fas', '.tree'))

            single_cmd_params = [msa_path, phylogenetic_raw_tree_path_without_bootstrap_values,
                                 pipeline_step_tmp_dir, pruned_tree_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=100,
                                                   job_name_suffix='adjust_tree', queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 25.	fix_msa.py
    # Input: (1) A path to an MSA file that (might) contain too few genomes and/or
    #                                                       ambiguous characters (i.e., other than a, c, g, t, -)
    #        (2) A path to a fixed msa
    # Output: a fixed msa at (2)
    step_number = '25'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_fixed_dna_msa'
    pipeline_step_output_dir_25, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        script_path = f'/bioseq/sincopa/fix_msa.py'
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for msa_file in os.listdir(pipeline_step_output_dir_23):
            # preparing parameters
            msa_path = os.path.join(pipeline_step_output_dir_23, msa_file)
            fixed_msa_output_path = os.path.join(pipeline_step_output_dir_25, msa_file)

            single_cmd_params = [msa_path, fixed_msa_output_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   num_of_cmds_per_job=250,
                                                   job_name_suffix='fix_msa', queue_name=args.queue_name)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 26.	compute_homoplasy.py
    # Input: (1) A path to an MSA file without ambiguous characters (only 'a', 'c', 'g', 't', and '-' are allowed)
    #        (2) A path to a corresponding background tree
    #        (3) A path to a control file for the homoplasy calculation cpp script
    #        (4) A path to a file in which the homoplasy will be written
    # Output: homoplasy of the input msa at (2)
    step_number = '26'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_homoplasy'
    pipeline_step_output_dir_26, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        script_path = f'/bioseq/sincopa/compute_homoplasy.py'
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for msa_file in os.listdir(pipeline_step_output_dir_23):
            # preparing parameters
            fixed_msa_path = os.path.join(pipeline_step_output_dir_25, msa_file)
            tree_path = os.path.join(pipeline_step_output_dir_24, msa_file.replace('fas', 'tree'))
            control_path = os.path.join(pipeline_step_tmp_dir, msa_file.replace('fas', 'control'))
            output_path = os.path.join(pipeline_step_output_dir_26, msa_file.replace('fas', 'homoplasy'))

            single_cmd_params = [fixed_msa_path, tree_path, control_path, output_path]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   job_name_suffix='homoplasy', queue_name=args.queue_name,
                                                   num_of_cmds_per_job=50,
                                                   required_modules_as_list=[consts.GCC])

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 27.	compute_sweeps_score.py
    # Input: (1) A path to an MSA file to be analyzed
    #        (2) A path to a homoplasy
    #        (3) A path to a file in which the sweep scores will be written
    #        (4) A path to a file in which the sweep scores plot will be saved
    # Output: sweep scores and sweep scores plot at (3) and (4), respectively
    step_number = '27'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_sweeps_scores_computation'
    pipeline_step_output_dir_27, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    window_size = 50
    if not os.path.exists(done_file_path):
        script_path = f'/bioseq/sincopa/compute_sweeps_score.py'
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for fixed_msa_file in os.listdir(pipeline_step_output_dir_25):
            # preparing parameters
            fixed_msa_path = os.path.join(pipeline_step_output_dir_25, fixed_msa_file)
            homoplasy_path = os.path.join(pipeline_step_output_dir_26, fixed_msa_file.replace('fas', 'homoplasy'))
            output_scores_path = os.path.join(pipeline_step_output_dir_27, fixed_msa_file.replace('fas', 'txt'))
            output_meta_path = os.path.join(pipeline_step_output_dir_27, fixed_msa_file.replace('fas', 'csv'))
            output_plot_path = os.path.join(pipeline_step_output_dir_27, fixed_msa_file.replace('fas', 'png'))

            single_cmd_params = [fixed_msa_path, homoplasy_path, output_scores_path,
                                 output_meta_path, output_plot_path, window_size]
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   job_name_suffix='sweeps', queue_name=args.queue_name,
                                                   num_of_cmds_per_job=100)

        wait_for_results(logger, times_logger, os.path.split(script_path)[-1], pipeline_step_tmp_dir,
                         num_of_batches, error_file_path, email=args.email)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # 28. analyse sweeps
    # Input: a path to a folder with sweeps score files
    # Output: concatenated csv file with all the analysis metadata
    # CANNOT be parallelized on cluster
    step_number = '28'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = 'sweeps_analysis'
    sweeps_analysis_dir = os.path.join(output_dir, f'{step_number}_{step_name}')
    os.makedirs(sweeps_analysis_dir, exist_ok=True)
    sweeps_summary_path = os.path.join(sweeps_analysis_dir, 'sweeps_summary.csv')
    sorted_sweeps_summary_path = os.path.join(sweeps_analysis_dir, 'sorted_sweeps_summary.csv')

    if not os.path.exists(sorted_sweeps_summary_path):
        logger.info('Concatenating sweeps analysis metadata results...')

        # write summary header
        header = 'msa_name,max_score,number_of_sequences,centrality,msa_length,window_size,index_of_max,' \
                 'mean_score,median_score,max_mean_division,max_median_division,relative_location_of_peak,' \
                 'apd,pi,above95,above75,above50,above25,above05'
        with open(sweeps_summary_path, 'w') as f:
            f.write(f'{header}\n')

        # add content
        execute(logger, f'cat {pipeline_step_output_dir_27}/*csv >> {sweeps_summary_path}', process_is_string=True)

        df = pd.read_csv(sweeps_summary_path)

        df.sort_values(by=['max_score', 'number_of_sequences', 'centrality', 'msa_length'], ascending=False,
                       inplace=True)
        df.to_csv(sorted_sweeps_summary_path, index=False)

        run_number2species = {'159716760987322741064374127401': 'C. trachomatis, r/m=0.66',
                              '159717709089300979385014858842': 'H. pylori, r/m=21.11',
                              '159717758554609942254414562338': 'S. pneumoniae, r/m=5.17',
                              '159717782526726941828676365733': 'N. meningitidis, r/m=14.62',
                              '159717817742964662930817004125': 'E. coli r/m=0.38',
                              '159717825610100634511315581997': 'M. tuberculosis, r/m=0'}

        for column in df.columns:
            if column == 'msa_name' or 'above' in column:
                logger.info(f'Skipping {column} column...')
                continue

            fig = plt.figure()
            plt.title(f'{column}\n({run_number2species.get(run_number, "unrecognized species")})')
            plt.hist(df[column], bins=50)  # TODO: density=True
            fig.savefig(f'{sweeps_analysis_dir}/{column}.png', dpi=500, bbox_inches='tight')
            plt.close()
            fig = plt.figure()
            plt.title(f'{column}\n({run_number2species.get(run_number, "unrecognized species")})')
            df[column].plot.kde()
            fig.savefig(f'{sweeps_analysis_dir}/{column}_smooth.png', dpi=500, bbox_inches='tight')
            plt.close()
            fig = plt.figure()
            percentile = 10 / 100
            plt.title(
                f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
            df[column][df[column] >= df[column].quantile(1 - percentile)].plot.kde()
            fig.savefig(f'{sweeps_analysis_dir}/{column}_smooth_top_{percentile * 100}_percentile.png', dpi=500,
                        bbox_inches='tight')
            plt.close()
            fig = plt.figure()
            percentile = 20 / 100
            plt.title(
                f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
            df[column][df[column] >= df[column].quantile(1 - percentile)].plot.kde()
            fig.savefig(f'{sweeps_analysis_dir}/{column}_smooth_top_{percentile * 100}_percentile.png', dpi=500,
                        bbox_inches='tight')
            plt.close()
            fig = plt.figure()
            percentile = 5 / 1000
            plt.title(
                f'{column} ({percentile} trimmed)\n({run_number2species.get(run_number, "unrecognized species")})')
            plt.hist(np.sort(df[column])[int(len(df) * percentile): len(df) - int(len(df) * percentile)], bins=50)
            fig.savefig(f'{sweeps_analysis_dir}/{column}_trimmed.png', dpi=500, bbox_inches='tight')
            plt.close()
            fig = plt.figure()
            percentile = 1 / 100
            plt.title(
                f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
            plt.hist(np.sort(df[column])[len(df) - int(len(df) * percentile):], bins=30)
            fig.savefig(f'{sweeps_analysis_dir}/{column}_top_{percentile * 100}_percentile.png', dpi=500,
                        bbox_inches='tight')
            plt.close()
            fig = plt.figure()
            percentile = 5 / 100
            plt.title(
                f'{column} (top {int(percentile * 100)}%)\n({run_number2species.get(run_number, "unrecognized species")})')
            plt.hist(np.sort(df[column])[len(df) - int(len(df) * percentile):], bins=30)
            fig.savefig(f'{sweeps_analysis_dir}/{column}_top_{percentile * 100}_percentile.png', dpi=500,
                        bbox_inches='tight')
            plt.close()

        # No need to wait...
    else:
        logger.info(f'done file {sorted_sweeps_summary_path} already exists. Skipping step...')

    # preventing folder's deletion (output dir is being deleted
    # logger.info(f'Moving {pipeline_step_output_dir_24} TO {os.path.join(meta_output_dir, step_name)}')
    # try:
    #     shutil.move(pipeline_step_output_dir_23, meta_output_dir)
    #     shutil.move(pipeline_step_output_dir_24, meta_output_dir)
    # except FileExistsError:
    #     pass


def cleanup_results(args, logger, meta_output_dir, final_output_dir_name, output_dir):
    logger.info('Cleaning up...')

    # remove intermediate results (including tmp_dir)
    remove_path(logger, output_dir)

    # remove raw data from the server
    for path_to_remove in [os.path.join(meta_output_dir, final_output_dir_name, x) for x in
                           ['12_orthologs_groups_dna_sequences',
                            '13_orthologs_groups_aa_sequences',
                            '14_aligned_aa_orthologs_groups',
                            '15_aligned_core_proteome',
                            '17_orfs_statistics',
                            '18_induce_dna_msa_by_aa_msa']]:
        # remove intermediate results
        remove_path(logger, path_to_remove)


def main(args):
    start_time = time()

    logger, times_logger, meta_output_dir, error_file_path, run_number, output_html_path, output_url, meta_output_url,\
        output_dir, tmp_dir, done_files_dir = prepare_pipeline_framework(args)

    try:
        validate_arguments(args)
        data_path, number_of_genomes = prepare_and_verify_input_data(args, logger, meta_output_dir, error_file_path, output_dir)

        orfs_dir, orfs_step_number, final_orthologs_table_file_path, phylogenetic_raw_tree_path, final_output_dir_name = \
            run_main_pipeline(args, logger, times_logger, meta_output_dir, error_file_path, run_number,
                              output_html_path, output_dir, tmp_dir, done_files_dir, data_path, number_of_genomes)

        # run_pipeline_extensions(args, logger, error_file_path, run_number, output_dir, tmp_dir, done_files_dir,
        #                         data_path,
        #                         orfs_dir, orfs_step_number, final_orthologs_table_file_path,
        #                         phylogenetic_raw_tree_path, final_output_dir_name)

        status = 'is done'

        if run_number.lower() != 'example' and 'oren' not in args.email and not consts.TEST:
            cleanup_results(args, logger, meta_output_dir, final_output_dir_name, output_dir)
    except Exception as e:
        status = 'was failed'
        report_error_in_main_pipeline_to_admin(logger, e, meta_output_dir, error_file_path, run_number, output_html_path,
                                               meta_output_url)

    total_time = int(time() - start_time)
    report_main_pipeline_result_to_user(args, logger, status, total_time, output_url, run_number)

    logger.info('Done.')


if __name__ == '__main__':
    print(f'sys.path is\n{sys.path}')
    arguments = get_arguments()

    main(arguments)
