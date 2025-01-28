import os
import sys
from sys import argv
import argparse
import traceback
import logging
import shutil
import subprocess
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger, get_job_times_logger, none_or_str, prepare_directories, \
    submit_batch, wait_for_results, add_results_to_final_dir
from auxiliaries.file_writer import write_to_file
from auxiliaries.logic_auxiliaries import define_intervals
from auxiliaries.infer_orthogroups_logic import infer_orthogroups
from auxiliaries import consts


def create_pseudo_genome_from_ogs(logger, times_logger, base_step_number, final_substep_number, error_file_path, output_dir,
                           tmp_dir, done_files_dir, final_orthogroups_file_path, subset_proteins_fasta_path,
                           max_parallel_jobs, genomes_batch_id, pseudo_genome_mode):
    final_orthogroups_df = pd.read_csv(final_orthogroups_file_path)
    number_of_ogs = len(final_orthogroups_df.index)

    step_number = f'{base_step_number}_{final_substep_number + 1}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_fasta'
    script_path = os.path.join(consts.SRC_DIR, 'steps/extract_orfs.py')
    orthogroups_fasta_dir_path, pipeline_step_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    orthologs_aa_dir_path = os.path.join(orthogroups_fasta_dir_path, 'orthogroups_aa')
    orthologs_aa_aligned_dir_path = os.path.join(orthogroups_fasta_dir_path, 'orthogroups_aa_aligned')
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')
        os.makedirs(orthologs_aa_dir_path, exist_ok=True)
        os.makedirs(orthologs_aa_aligned_dir_path, exist_ok=True)

        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs

        ogs_intervals = define_intervals(0, number_of_ogs - 1, max_parallel_jobs)
        for og_number_start, og_number_end in ogs_intervals:
            single_cmd_params = [None,
                                 subset_proteins_fasta_path,
                                 final_orthogroups_file_path,
                                 og_number_start,
                                 og_number_end,
                                 None,
                                 orthologs_aa_dir_path,
                                 orthologs_aa_aligned_dir_path if pseudo_genome_mode == 'consensus_gene' else None,
                                 None]
            all_cmds_params.append(single_cmd_params)

        num_of_batches = submit_batch(logger, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                                   error_file_path,
                                                   num_of_cmds_per_job=1,
                                                   job_name_suffix='orfs_extraction',
                                                   queue_name=args.queue_name,
                                                   account_name=args.account_name,
                                                   node_name=args.node_name)

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, num_of_batches, error_file_path)
        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    if pseudo_genome_mode == 'consensus_gene':
        step_number = f'{base_step_number}_{final_substep_number + 2}'
        logger.info(f'Step {step_number}: {"_" * 100}')
        step_name = f'{step_number}_ogs_consensus'
        script_path = os.path.join(consts.SRC_DIR, 'steps/calculate_consensus_from_msa.py')
        ogs_consensus_dir_path, ogs_consensus_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
        done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
        if not os.path.exists(done_file_path):
            logger.info('Calculating consensus sequences for all aligned OGs...')

            job_input_files = os.path.join(ogs_consensus_tmp_dir, 'job_input_files')
            os.makedirs(job_input_files, exist_ok=True)

            job_index_to_ogs_paths = defaultdict(list)
            for i, og_aligned_file_name in enumerate(os.listdir(orthologs_aa_aligned_dir_path)):
                og_aligned_path = os.path.join(orthologs_aa_aligned_dir_path, og_aligned_file_name)
                job_index = i % max_parallel_jobs
                job_index_to_ogs_paths[job_index].append(og_aligned_path)

            all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
            for job_index, ogs_paths in job_index_to_ogs_paths.items():
                job_input_path = os.path.join(job_input_files, f'job_input_{job_index}.txt')
                with open(job_input_path, 'w') as f:
                    f.write('\n'.join(ogs_paths))

                single_cmd_params = [job_input_path, ogs_consensus_dir_path]
                all_cmds_params.append(single_cmd_params)

            num_of_batches = submit_batch(logger, script_path, all_cmds_params, ogs_consensus_tmp_dir, error_file_path,
                                          num_of_cmds_per_job=1, job_name_suffix='ogs_consensus',
                                          queue_name=args.queue_name, account_name=args.account_name,
                                          node_name=args.node_name)

            wait_for_results(logger, times_logger, step_name, ogs_consensus_tmp_dir, num_of_batches, error_file_path)
            write_to_file(logger, done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists. Skipping step...')


    step_number = f'{base_step_number}_{final_substep_number + 3}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_pseudo_genome'
    pseudo_genome_dir_path, pseudo_genome_tmp_dir = prepare_directories(logger, output_dir, tmp_dir, step_name)
    done_file_path = os.path.join(done_files_dir, f'{step_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Creating pseudo genome for batch...')

        records = []
        record_to_og = {}

        if pseudo_genome_mode == 'first_gene':
            for filename in os.listdir(orthologs_aa_dir_path):
                if not filename.endswith('.faa'):
                    continue

                og_name = filename.split('.')[0]
                og_path = os.path.join(orthologs_aa_dir_path, filename)

                # read only the first sequence from each og
                first_record = SeqIO.parse(og_path, 'fasta').__next__()
                first_record.id = f'pseudo_genome_{genomes_batch_id}:{og_name}_representative'

                records.append(first_record)
                record_to_og[first_record.id] = og_name

        elif pseudo_genome_mode == 'consensus_gene':
            for filename in os.listdir(ogs_consensus_dir_path):
                if not filename.endswith('.faa'):
                    continue

                og_name = filename.split('.')[0]
                og_path = os.path.join(ogs_consensus_dir_path, filename)

                record = SeqIO.parse(og_path, 'fasta').__next__()
                record.id = f'pseudo_genome_{genomes_batch_id}:{og_name}_consensus'

                records.append(record)
                record_to_og[record.id] = og_name

        pseudo_genome_fasta_path = os.path.join(pseudo_genome_dir_path, f'pseudo_genome_{genomes_batch_id}.faa')
        SeqIO.write(records, pseudo_genome_fasta_path, 'fasta')
        logger.info(f'Wrote {len(records)} records to {pseudo_genome_fasta_path} (pseudo genome type is "{pseudo_genome_mode}")')

        record_to_og_df = pd.DataFrame(record_to_og.items(), columns=['representative_gene', 'OG_name'])
        orthogroups_with_representative_df = final_orthogroups_df.merge(record_to_og_df, on='OG_name')

        orthogroups_with_representative_path = os.path.join(pseudo_genome_dir_path, f'orthogroups_with_representative_{genomes_batch_id}.csv')
        orthogroups_with_representative_df.to_csv(orthogroups_with_representative_path, index=False)
        logger.info(f'Wrote orthogroups with representative to {orthogroups_with_representative_path}')

        write_to_file(logger, done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def infer_orthogroups_on_genomes_batch(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir,
            done_files_dir, translated_orfs_dir, genomes_names_path,
            queue_name, account_name, node_name, identity_cutoff, coverage_cutoff,
            e_value_cutoff, sensitivity, max_parallel_jobs, genomes_batch_id, pseudo_genome_mode, run_optimized_mmseqs, use_parquet,
            verbose, add_orphan_genes_to_ogs):

    with open(genomes_names_path, 'r') as genomes_names_fp:
        genomes_names = genomes_names_fp.read().split('\n')

    subset_proteomes_dir = os.path.join(output_dir, 'proteomes')
    if not os.path.exists(subset_proteomes_dir):
        os.makedirs(subset_proteomes_dir, exist_ok=True)
        for genome_name in genomes_names:
            shutil.copy(os.path.join(translated_orfs_dir, f'{genome_name}.faa'),
                        os.path.join(subset_proteomes_dir, f'{genome_name}.faa'))

    subset_proteins_fasta_path = os.path.join(output_dir, 'all_proteomes.faa')
    if not os.path.exists(subset_proteins_fasta_path):
        cmd = f"cat {os.path.join(subset_proteomes_dir, '*')} > {subset_proteins_fasta_path}"
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True)

    final_orthogroups_dir_path, orphan_genes_dir, final_substep_number = infer_orthogroups(
        logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir,
            done_files_dir, subset_proteomes_dir, subset_proteins_fasta_path, genomes_names_path,
            queue_name, account_name, node_name, identity_cutoff, coverage_cutoff,
            e_value_cutoff, sensitivity, max_parallel_jobs, run_optimized_mmseqs, use_parquet,
            verbose, add_orphan_genes_to_ogs)

    final_orthogroups_file_path = os.path.join(final_orthogroups_dir_path, 'orthogroups.csv')

    create_pseudo_genome_from_ogs(logger, times_logger, base_step_number, final_substep_number, error_file_path, output_dir,
                           tmp_dir, done_files_dir, final_orthogroups_file_path, subset_proteins_fasta_path,
                           max_parallel_jobs, genomes_batch_id, pseudo_genome_mode)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('step_number', help='')
    parser.add_argument('output_dir', help='')
    parser.add_argument('tmp_dir', help='')
    parser.add_argument('done_files_dir', help='')
    parser.add_argument('translated_orfs_dir', help='')
    parser.add_argument('genomes_names_path', help='')
    parser.add_argument('queue_name', help='')
    parser.add_argument('account_name', help='')
    parser.add_argument('node_name', type=none_or_str, help='')
    parser.add_argument('identity_cutoff', help='', type=float)
    parser.add_argument('coverage_cutoff', help='', type=float)
    parser.add_argument('e_value_cutoff', help='', type=float)
    parser.add_argument('sensitivity', type=float)
    parser.add_argument('max_parallel_jobs', help='', type=int)
    parser.add_argument('genomes_batch_id', help='', type=int)
    parser.add_argument('pseudo_genome_mode', help='', type=str)
    parser.add_argument('--run_optimized_mmseqs', help='', action='store_true')
    parser.add_argument('--use_parquet', action='store_true')
    parser.add_argument('--add_orphan_genes_to_ogs', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)
    times_logger = get_job_times_logger(args.logs_dir)

    logger.info(script_run_message)
    try:
        infer_orthogroups_on_genomes_batch(
            logger, times_logger, args.step_number, args.error_file_path, args.output_dir, args.tmp_dir,
            args.done_files_dir, args.translated_orfs_dir, args.genomes_names_path,
            args.queue_name, args.account_name, args.node_name, args.identity_cutoff, args.coverage_cutoff,
            args.e_value_cutoff, args.sensitivity, args.max_parallel_jobs, args.genomes_batch_id,
            args.pseudo_genome_mode, args.run_optimized_mmseqs, args.use_parquet, args.verbose,
            args.add_orphan_genes_to_ogs)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
