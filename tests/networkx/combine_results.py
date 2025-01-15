import pandas as pd
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    outputs_dir = os.path.join(SCRIPT_DIR, 'outputs')
    dfs = []

    for dir_name in os.listdir(outputs_dir):
        if not os.path.isdir(os.path.join(outputs_dir, dir_name)):
            continue
        summary_path = os.path.join(outputs_dir, dir_name, 'summary.csv')
        df = pd.read_csv(summary_path)
        dfs.append(df)

    combined_df = pd.concat(dfs)
    combined_df.to_csv(os.path.join(outputs_dir, 'combined_summary.csv'), index=False)


if __name__ == '__main__':
    main()
