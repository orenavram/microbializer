import os
import sys
from sys import argv
import argparse
import traceback
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger, prepare_directories, submit_batch, \
    wait_for_results, submit_mini_batch, get_job_times_logger
from auxiliaries.cluster_mmseqs_hits_to_orthogroups import cluster_mmseqs_hits_to_orthogroups, run_mmseqs_and_extract_hits


def full_orthogroups_infernece(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir, done_files_dir,
                               translated_orfs_dir, all_proteins_path, strains_names_path, queue_name,
                               account_name, identity_cutoff, coverage_cutoff, e_value_cutoff, max_parallel_jobs,
                               run_optimized_mmseqs, unify_clusters_after_mmseqs, use_parquet, prepare_mcl_v2,
                               run_mcl_on_all_hits_together, use_linux_to_parse_big_files, verbose):

    orthologs_output_dir, paralogs_output_dir, orthologs_scores_statistics_dir, paralogs_scores_statistics_dir = \
        run_mmseqs_and_extract_hits(logger, times_logger, base_step_number, error_file_path, output_dir, tmp_dir,
                                    done_files_dir, translated_orfs_dir, all_proteins_path, strains_names_path, queue_name,
                                    account_name, identity_cutoff, coverage_cutoff, e_value_cutoff, max_parallel_jobs,
                                    run_optimized_mmseqs, use_parquet, use_linux_to_parse_big_files, verbose)

    if unify_clusters_after_mmseqs:
        return

    cluster_mmseqs_hits_to_orthogroups(logger, times_logger, error_file_path, output_dir, tmp_dir, done_files_dir,
                                       orthologs_output_dir, orthologs_scores_statistics_dir, paralogs_output_dir,
                                       paralogs_scores_statistics_dir, max_parallel_jobs, base_step_number,
                                       4, account_name, queue_name, use_parquet, prepare_mcl_v2, strains_names_path,
                                       run_mcl_on_all_hits_together)


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
    parser.add_argument('max_parallel_jobs', help='', type=int)
    parser.add_argument('--run_optimized_mmseqs', help='', action='store_true')
    parser.add_argument('--unify_clusters_after_mmseqs', help='', action='store_true')
    parser.add_argument('--use_parquet', action='store_true')
    parser.add_argument('--prepare_mcl_v2', action='store_true')
    parser.add_argument('--run_mcl_on_all_hits_together', action='store_true')
    parser.add_argument('--use_linux_to_parse_big_files', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)
    times_logger = get_job_times_logger(args.logs_dir)

    logger.info(script_run_message)
    try:
        full_orthogroups_infernece(logger, times_logger, args.step_number, args.error_file_path, args.output_dir,
                                   args.tmp_dir, args.done_files_dir, args.translated_orfs_dir, args.all_proteins_path,
                                   args.strains_names_path, args.queue_name, args.account_name,
                                   args.identity_cutoff, args.coverage_cutoff, args.e_value_cutoff,
                                   args.max_parallel_jobs, args.run_optimized_mmseqs, args.unify_clusters_after_mmseqs,
                                   args.use_parquet, args.prepare_mcl_v2, args.run_mcl_on_all_hits_together,
                                   args.use_linux_to_parse_big_files, args.verbose)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
