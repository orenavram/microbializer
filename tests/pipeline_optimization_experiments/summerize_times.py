import pandas as pd
from pathlib import Path
import re

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / '73_ecoli_optimizations_summary'

RELEVANT_LOG_STEPS = {
    ('no', 'no'): ['05_1_all_vs_all_analysis', '05_2_max_rbh_scores', '05_3_mmseqs_paralogs', '05_4_normalize_scores',
                   '05_5_concatenate_hits', '05_6_putative_table', '05_7_mcl_input_files', '05_8_mcl_analysis',
                   '05_9_verified_clusters', '05_10_verified_table'],
    ('no', 'yes'): ['05_1_all_vs_all_analysis', '05_2_extract_rbh_hits', '05_3_paralogs', '05_4_normalize_scores',
                    '05_5_concatenate_hits', '05_6_putative_table', '05_7_mcl_input_files', '05_8_mcl_analysis',
                    '05_9_verified_clusters', '05_10_verified_table'],
    ('mergeAtEnd', 'no'): ['04_1_cluster_proteomes', '04_2_filter_fastas_by_cluster', '05_infer_orthogroups'],
    ('mergeAtEnd', 'yes'): ['04_1_cluster_proteomes', '04_2_filter_fastas_by_cluster', '05_infer_orthogroups'],
    ('mergeAfterMmseqs', 'no'): ['04_1_cluster_proteomes', '04_2_filter_fastas_by_cluster', '05_1_infer_orthogroups',
                                 '05_2_concat_clusters_rbhs', '05_3_concat_clusters_paralogs', '05_4_normalize_scores',
                                 '05_5_concatenate_hits', '05_6_putative_table', '05_7_mcl_input_files',
                                 '05_8_mcl_analysis', '05_9_verified_clusters', '05_10_verified_table'],
    ('mergeAfterMmseqs', 'yes'): ['04_1_cluster_proteomes', '04_2_filter_fastas_by_cluster', '05_1_infer_orthogroups',
                                  '05_2_concat_clusters_rbhs', '05_3_concat_clusters_paralogs', '05_4_normalize_scores',
                                  '05_5_concatenate_hits', '05_6_putative_table', '05_7_mcl_input_files',
                                  '05_8_mcl_analysis', '05_9_verified_clusters', '05_10_verified_table']
}


def main():
    experiments_times = {}
    for dir in DATA_DIR.iterdir():
        if not dir.is_dir():
            continue

        dir_name_match = re.search(r'outputs_preCluster_(.*)_optimizedMmseqs_(.*)_useParquet_(.*)', dir.name)
        preCluster = dir_name_match.group(1)
        optimizedMmseqs = dir_name_match.group(2)
        useParquet = dir_name_match.group(3)

        relevant_log_steps = RELEVANT_LOG_STEPS[(preCluster, optimizedMmseqs)]
        times_log_path = dir / 'times_log.txt'
        if not times_log_path.exists():
            experiments_times[(preCluster, optimizedMmseqs, useParquet)] = pd.NA
            continue

        with open(times_log_path) as f:
            times_file_lines = f.readlines()

        times = []
        for line in times_file_lines:
            if any(step in line for step in relevant_log_steps):
                time_match = re.search(r'cpus used per job = (.*) wallclock time', line)
                if time_match:
                    time_string = time_match.group(1)
                else:
                    time_match = re.search(r'cumulatively they took (.*) wallclock time', line)
                    time_string = time_match.group(1)

                times.append(pd.Timedelta(time_string))

        total_time = sum(times, pd.Timedelta(0))
        experiments_times[(preCluster, optimizedMmseqs, useParquet)] = total_time

    experiments_times_df = pd.DataFrame(
        [(key[0], key[1], key[2], value) for key, value in experiments_times.items()],
        columns=['preCluster', 'optimizedMmseqs', 'useParquet', 'totalTime']
    )
    experiments_times_df.to_csv(DATA_DIR / 'experiments_times.csv', index=False)


if __name__ == '__main__':
    main()
