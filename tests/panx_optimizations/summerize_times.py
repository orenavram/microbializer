import pandas as pd
from pathlib import Path
import re

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / 'benchmark_1060_c_001'

RELEVANT_STEPS_PREFIXES = {
    'original': ['05_'],
    'pseudo_genomes': ['05_subsets_inference', '06_'],
    'pseudo_genomes_consensus': ['05_subsets_inference', '06_'],
}


def main():
    experiments_times = {}
    for dir in DATA_DIR.iterdir():
        if not dir.is_dir():
            continue

        relevant_log_steps = RELEVANT_STEPS_PREFIXES[dir.name]
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
        experiments_times[dir.name] = total_time

    experiments_times_df = pd.DataFrame(list(experiments_times.items()), columns=['inference_method', 'totalTime'])
    experiments_times_df.to_csv(DATA_DIR / 'experiments_times.csv', index=False)


if __name__ == '__main__':
    main()
