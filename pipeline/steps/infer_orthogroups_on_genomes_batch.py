import sys
import argparse
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_times_logger, prepare_directories, run_step, \
    submit_batch, wait_for_results, add_results_to_final_dir, add_default_step_args, write_done_file, str_to_bool
from auxiliaries.logic_auxiliaries import split_ogs_to_jobs_inputs_files_by_og_sizes
from auxiliaries.infer_orthogroups_logic import infer_orthogroups
from auxiliaries.configuration import Config
from auxiliaries import consts


def create_pseudo_genome_from_ogs(
        logger, times_logger, config, base_step_number, previous_substep_number, final_orthogroups_file_path,
        subset_proteins_fasta_path, genomes_batch_id):
    final_orthogroups_df = pd.read_csv(final_orthogroups_file_path)

    step_number = f'{base_step_number}_{previous_substep_number + 1}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_orthogroups_fasta'
    script_path = consts.SRC_DIR / 'steps' / 'extract_orfs.py'
    orthogroups_fasta_dir_path, pipeline_step_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, step_name)
    orthologs_aa_dir_path = orthogroups_fasta_dir_path / 'orthogroups_aa'
    orthologs_aa_aligned_dir_path = orthogroups_fasta_dir_path / 'orthogroups_aa_aligned'
    orthologs_aa_consensus_dir_path = orthogroups_fasta_dir_path / 'orthogroups_aa_consensus'
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')
        orthologs_aa_dir_path.mkdir(parents=True, exist_ok=True)
        orthologs_aa_aligned_dir_path.mkdir(parents=True, exist_ok=True)
        orthologs_aa_consensus_dir_path.mkdir(parents=True, exist_ok=True)

        job_paths = split_ogs_to_jobs_inputs_files_by_og_sizes(
            final_orthogroups_df, pipeline_step_tmp_dir, config.max_parallel_jobs)
        all_cmds_params = []
        for job_path in job_paths:
            single_cmd_params = [None,
                                 subset_proteins_fasta_path,
                                 final_orthogroups_file_path,
                                 job_path,
                                 None,
                                 orthologs_aa_dir_path,
                                 orthologs_aa_aligned_dir_path if config.pseudo_genome_mode == 'consensus_gene' else None,
                                 None,
                                 orthologs_aa_consensus_dir_path if config.pseudo_genome_mode == 'consensus_gene' else None]
            all_cmds_params.append(single_cmd_params)

        num_of_batches = submit_batch(logger, config, script_path, all_cmds_params, pipeline_step_tmp_dir,
                                      'orfs_extraction')

        wait_for_results(logger, times_logger, step_name, pipeline_step_tmp_dir, num_of_batches, config.error_file_path)
        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')

    step_number = f'{base_step_number}_{previous_substep_number + 2}'
    logger.info(f'Step {step_number}: {"_" * 100}')
    step_name = f'{step_number}_pseudo_genome'
    pseudo_genome_dir_path, pseudo_genome_tmp_dir = prepare_directories(
        logger, config.steps_results_dir, config.tmp_dir, step_name)
    done_file_path = config.done_files_dir / f'{step_name}.txt'
    if not done_file_path.exists():
        logger.info('Creating pseudo genome for batch...')

        records = []
        record_to_og = {}

        if config.pseudo_genome_mode == 'first_gene':
            for og_path in orthologs_aa_dir_path.glob('*.faa'):
                og_name = og_path.stem

                # read only the first sequence from each og
                first_record = SeqIO.parse(og_path, 'fasta').__next__()
                first_record.id = f'pseudo_genome_{genomes_batch_id}:{og_name}_representative'

                records.append(first_record)
                record_to_og[first_record.id] = og_name

        elif config.pseudo_genome_mode == 'consensus_gene':
            for og_path in orthologs_aa_consensus_dir_path.glob('*.faa'):
                og_name = og_path.stem

                record = SeqIO.parse(og_path, 'fasta').__next__()
                record.id = f'pseudo_genome_{genomes_batch_id}:{og_name}_consensus'

                records.append(record)
                record_to_og[record.id] = og_name

        pseudo_genome_fasta_path = pseudo_genome_dir_path / f'pseudo_genome_{genomes_batch_id}.faa'
        SeqIO.write(records, pseudo_genome_fasta_path, 'fasta')
        logger.info(f'Wrote {len(records)} records to {pseudo_genome_fasta_path} (pseudo genome type is '
                    f'"{config.pseudo_genome_mode}")')

        record_to_og_df = pd.DataFrame(record_to_og.items(), columns=['representative_gene', 'OG_name'])
        orthogroups_with_representative_df = final_orthogroups_df.merge(record_to_og_df, on='OG_name')

        orthogroups_with_representative_path = pseudo_genome_dir_path / f'orthogroups_with_representative_{genomes_batch_id}.csv'
        orthogroups_with_representative_df.to_csv(orthogroups_with_representative_path, index=False)
        logger.info(f'Wrote orthogroups with representative to {orthogroups_with_representative_path}')

        write_done_file(logger, done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists. Skipping step...')


def infer_orthogroups_on_genomes_batch(logger, times_logger, config, step_number, translated_orfs_dir,
                                       genomes_batch_id):
    with open(config.genomes_names_path, 'r') as genomes_names_fp:
        genomes_names = genomes_names_fp.read().split('\n')

    subset_proteomes_dir = config.steps_results_dir / 'proteomes'
    if not subset_proteomes_dir.exists():
        subset_proteomes_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f'Copying translated ORFs to {subset_proteomes_dir}...')
        for genome_name in genomes_names:
            shutil.copy(translated_orfs_dir / f'{genome_name}.faa', subset_proteomes_dir / f'{genome_name}.faa')
            logger.info(f'Copied {genome_name}.faa to {subset_proteomes_dir}')

    subset_proteins_fasta_path = config.steps_results_dir / 'all_proteomes.faa'
    if not subset_proteins_fasta_path.exists():
        cmd = f"cat {subset_proteomes_dir / '*'} > {subset_proteins_fasta_path}"
        logger.info(f'Calling: {cmd}')
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

    final_orthogroups_dir_path, orphan_genes_dir, final_substep_number = infer_orthogroups(
        logger, times_logger, config, step_number, subset_proteomes_dir, subset_proteins_fasta_path)

    final_orthogroups_file_path = final_orthogroups_dir_path / 'orthogroups.csv'

    create_pseudo_genome_from_ogs(
        logger, times_logger, config, step_number, final_substep_number, final_orthogroups_file_path,
        subset_proteins_fasta_path, genomes_batch_id)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_path', help='')
    parser.add_argument('step_number', help='')
    parser.add_argument('translated_orfs_dir', type=Path, help='')
    parser.add_argument('genomes_batch_id', help='', type=int)
    add_default_step_args(parser)
    args = parser.parse_args()

    times_logger = get_job_times_logger(args.logs_dir, args.job_name, args.verbose)
    config = Config.from_csv(args.config_path)

    run_step(args, infer_orthogroups_on_genomes_batch, times_logger, config, args.step_number, args.translated_orfs_dir,
             args.genomes_batch_id)
