import os
import sys
import itertools
from sys import argv
import argparse
import traceback
from collections import defaultdict
import dask.dataframe as dd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, prepare_directories, submit_batch, \
    wait_for_results, submit_mini_batch, get_job_times_logger
from auxiliaries.logic_auxiliaries import define_intervals, add_score_column_to_mmseqs_output
from auxiliaries import consts
from auxiliaries.file_writer import write_to_file


def run_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                       all_proteins_path, strains_names, queue_name, account_name, identity_cutoff, coverage_cutoff,
                       e_value_cutoff, n_jobs_per_step):
    # a.	mmseqs2_all_vs_all.py
    step_number = f'{base_step_number}a'
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
                  f'--cpus {consts.MMSEQS_NUM_OF_CORES}']

        submit_mini_batch(logger, script_path, [params], pipeline_step_tmp_dir, error_file_path, queue_name,
                          account_name, job_name='mmseqs', num_of_cpus=consts.MMSEQS_NUM_OF_CORES,
                          memory=consts.MMSEQS_REQUIRED_MEMORY_GB, time_in_hours=72)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, 1, error_file_path)

        logger.info(f"Starting to read {m8_raw_output_path} and create a reduced version of it...")
        m8_df = dd.read_csv(m8_raw_output_path, sep='\t', names=consts.MMSEQS_OUTPUT_HEADER, dtype=consts.MMSEQS_OUTPUT_COLUMNS_TYPES)
        m8_df = m8_df[m8_df['query'] != m8_df['target']]
        add_score_column_to_mmseqs_output(m8_df)
        m8_df['query_genome'] = m8_df['query'].str.split(':').str[0]
        m8_df['target_genome'] = m8_df['target'].str.split(':').str[0]
        m8_df = m8_df[['query', 'query_genome', 'target', 'target_genome', 'score']]
        m8_df.to_parquet(m8_output_path)
        logger.info(f"{m8_output_path} was created successfully.")

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # b.	extract_rbh_hits.py
    step_number = f'{base_step_number}b'
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
                all_cmds_params.append(params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                       error_file_path,
                                                       num_of_cmds_per_job=1, # = max(1, len(all_cmds_params) // n_jobs_per_step),
                                                       job_name_suffix='rbh_analysis',
                                                       queue_name=queue_name,
                                                       account_name=account_name,
                                                       memory=consts.MMSEQS_PARSING_REQUIRED_MEMORY_GB)

            wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                             num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    # c. extract_paralogs.py
    step_number = f'{base_step_number}c'
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
            all_cmds_params.append(single_cmd_params)

        num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   error_file_path,
                                                   num_of_cmds_per_job=1, # = max(1, len(all_cmds_params) // n_jobs_per_step),
                                                   job_name_suffix='paralogs_analysis',
                                                   queue_name=queue_name,
                                                   account_name=account_name,
                                                   memory=consts.MMSEQS_PARSING_REQUIRED_MEMORY_GB)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir,
                         num_of_batches, error_file_path)

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    return orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir


def run_non_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                           translated_orfs_dir, strains_names, queue_name, account_name, identity_cutoff,
                           coverage_cutoff, e_value_cutoff, n_jobs_per_step):
    # a.	mmseqs2_all_vs_all.py
    # Input: (1) 2 input paths for 2 (different) genome files (query and target), g1 and g2
    #        (2) an output file path (with a suffix as follows: i_vs_j.tsv. especially relevant for the wrapper).
    # Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
    # Precisely, for each gene x in g1, query x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
    # Can be parallelized on cluster
    step_number = f'{base_step_number}a'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_all_vs_all_analysis'
    script_path = os.path.join(consts.SRC_DIR, 'steps/mmseqs_v1/mmseqs2_all_vs_all.py')
    orthologs_output_dir, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orthologs_scores_statistics_dir = os.path.join(orthologs_output_dir, 'scores_statistics')
    os.makedirs(orthologs_scores_statistics_dir, exist_ok=True)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        orfs_files = [filename for filename in os.listdir(translated_orfs_dir) if
                      filename.endswith('02d_translated_orfs')]
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
                all_cmds_params.append(single_cmd_params)

            num_of_batches, example_cmd = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                       error_file_path,
                                                       num_of_cmds_per_job=max(1,
                                                                               len(all_cmds_params) // n_jobs_per_step),
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

    # c. mmseqs2_paralogs.py
    # Input: max orthologs score for each gene
    # Output: paralogs in each genome
    step_number = f'{base_step_number}c'
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


def full_orthogroups_infernece(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                               translated_orfs_dir, all_proteins_path, strains_names_path, queue_name,
                               account_name, identity_cutoff, coverage_cutoff, e_value_cutoff, num_of_times_script_called,
                               run_optimized_mmseqs):
    with open(strains_names_path) as f:
        strains_names = f.read().rstrip().split('\n')
    n_jobs_per_step = consts.MAX_PARALLEL_JOBS // num_of_times_script_called

    if run_optimized_mmseqs:
        run_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir,
                           done_files_dir, all_proteins_path, strains_names, queue_name, account_name,
                           identity_cutoff, coverage_cutoff,  e_value_cutoff, n_jobs_per_step)
    else:
        run_non_unified_mmseqs(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir,
                               done_files_dir, translated_orfs_dir, strains_names, queue_name, account_name,
                               identity_cutoff, coverage_cutoff, e_value_cutoff, n_jobs_per_step)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('step_number', help='')
    parser.add_argument('output_dir', help='')
    parser.add_argument('tmp_dir', help='')
    parser.add_argument('done_files_dir', help='')
    parser.add_argument('translated_orfs_dir', help='')
    parser.add_argument('all_proteins_path', help='')
    parser.add_argument('strains_names_path', help='')
    parser.add_argument('queue_name', help='')
    parser.add_argument('account_name', help='')
    parser.add_argument('identity_cutoff', help='', type=float)
    parser.add_argument('coverage_cutoff', help='', type=float)
    parser.add_argument('e_value_cutoff', help='', type=float)
    parser.add_argument('num_of_times_script_called', help='', type=int)
    parser.add_argument('--run_optimized_mmseqs', help='', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    logger = get_job_logger(args.logs_dir)
    times_logger = get_job_times_logger(args.logs_dir)

    logger.info(script_run_message)
    try:
        full_orthogroups_infernece(logger, times_logger, args.step_number, args.error_file_path, args.output_dir,
                                   args.tmp_dir, args.done_files_dir, args.translated_orfs_dir, args.all_proteins_path,
                                   args.strains_names_path, args.queue_name, args.account_name,
                                   args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff,
                                   args.num_of_times_script_called, args.run_optimized_mmseqs)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
