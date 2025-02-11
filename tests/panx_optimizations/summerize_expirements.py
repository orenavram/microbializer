import pandas as pd
from pathlib import Path
import re

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / 'salmonella_300_c_001'

RELEVANT_STEPS_PREFIXES = {
    'original': ['05_'],
    'pseudo': ['05_subsets_inference', '06_'],
}


def main():
    experiments_info = []
    for dir in DATA_DIR.iterdir():
        if not dir.is_dir():
            continue

        if dir.name == 'original':
            relevant_log_steps = RELEVANT_STEPS_PREFIXES['original']
        else:
            relevant_log_steps = RELEVANT_STEPS_PREFIXES['pseudo']
        times_log_path = dir / 'times_log.txt'

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
                    if time_match:
                        time_string = time_match.group(1)
                    else:
                        time_match = re.search(r'took (.*).', line)
                        time_string = time_match.group(1)

                times.append(pd.Timedelta(time_string))

        total_time = sum(times, pd.Timedelta(0))

        orthogroups_path = dir / 'orthogroups.csv'
        orthogroups_df = pd.read_csv(orthogroups_path, dtype=str)
        num_orthogroups = len(orthogroups_df)

        scores_path = dir.glob('comparison_scores*.csv').__next__()
        scores_df = pd.read_csv(scores_path)
        scores = list(scores_df.iloc[0])

        experiments_info.append((dir.name, total_time, num_orthogroups, *scores))

    experiments_times_df = pd.DataFrame(experiments_info, columns=['inference_method', 'totalTime', 'num_orthogroups',
                                                                   'adjusted_rand_score', 'homogeneity', 'completeness',
                                                                   'v_measure', 'fowlkes_mallows_score'])
    experiments_times_df.to_csv(DATA_DIR / 'experiments_results.csv', index=False)


if __name__ == '__main__':
    main()
