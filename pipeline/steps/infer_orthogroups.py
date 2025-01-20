import os
import sys
from sys import argv
import argparse
import traceback
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger, get_job_times_logger, none_or_str
from auxiliaries.infer_orthogroups_logic import infer_orthogroups


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
    parser.add_argument('node_name', type=none_or_str, help='')
    parser.add_argument('identity_cutoff', help='', type=float)
    parser.add_argument('coverage_cutoff', help='', type=float)
    parser.add_argument('e_value_cutoff', help='', type=float)
    parser.add_argument('max_parallel_jobs', help='', type=int)
    parser.add_argument('--run_optimized_mmseqs', help='', action='store_true')
    parser.add_argument('--use_parquet', action='store_true')
    parser.add_argument('--use_linux_to_parse_big_files', action='store_true')
    parser.add_argument('--mmseqs_use_dbs', action='store_true')
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
        infer_orthogroups(
            logger, times_logger, args.step_number, args.error_file_path, args.output_dir, args.tmp_dir,
            args.done_files_dir, args.translated_orfs_dir, args.all_proteins_path, args.strains_names_path,
            args.queue_name, args.account_name, args.node_name, args.identity_cutoff, args.coverage_cutoff,
            args.e_value_cutoff, args.max_parallel_jobs, args.run_optimized_mmseqs, args.use_parquet,
            args.use_linux_to_parse_big_files, args.mmseqs_use_dbs, args.verbose, args.add_orphan_genes_to_ogs)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
