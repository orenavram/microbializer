import os
import sys
from sys import argv
import argparse
import logging
import subprocess
import traceback
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import fail, get_job_logger


def concatenate_hits(logger, clusters_dir, sub_dir_name, hits_filename, output_dir):
    hits_file_prefix, _ = os.path.splitext(hits_filename)
    output_hits_file = os.path.join(output_dir, hits_filename)
    output_statistics_file = os.path.join(output_dir, 'scores_statistics', f'{hits_file_prefix}.stats')
    all_statistics = []
    first_hits_file = True
    for cluster_dir in os.listdir(clusters_dir):
        hits_file_path = os.path.join(clusters_dir, cluster_dir, sub_dir_name, hits_filename)
        if os.path.exists(hits_file_path):
            if first_hits_file:
                cmd = f"cat {hits_file_path} >> {output_hits_file}"
                first_hits_file = False
            else:
                cmd = f"tail -n +2 {hits_file_path} >> {output_hits_file}"
            logger.info(f'Calling:\n{cmd}')
            subprocess.run(cmd, shell=True)

        statistics_file_path = os.path.join(clusters_dir, cluster_dir, sub_dir_name, 'scores_statistics', f'{hits_file_prefix}.stats')
        if os.path.exists(statistics_file_path):
            with open(statistics_file_path, 'r') as statistics_file:
                statistics = json.load(statistics_file)
            all_statistics.append(statistics)

    total_sum = 0
    total_records = 0
    for statistics_obj in all_statistics:
        total_sum += statistics_obj['sum']
        total_records += statistics_obj['number of records']

    overall_mean = total_sum / total_records if total_records > 0 else 0
    combined_json = {'mean': overall_mean, 'sum': total_sum, 'number of records': total_records}
    with open(output_statistics_file, 'w') as output_statistics_fp:
        json.dump(combined_json, output_statistics_fp)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('clusters_dir', help='dir of all clusters hits')
    parser.add_argument('sub_dir_name', help='name of relevant dir in the cluster dir')
    parser.add_argument('hits_filename', help='name of the hits file name')
    parser.add_argument('output_dir', help='output dir of concatenated files')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    parser.add_argument('--error_file_path', help='path to error file')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)
    try:
        concatenate_hits(logger, args.clusters_dir, args.sub_dir_name, args.hits_filename, args.output_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
        with open(args.error_file_path, 'a+') as f:
            traceback.print_exc(file=f)
