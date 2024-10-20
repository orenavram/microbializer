import subprocess
import os
import argparse
import pandas as pd


def concatenate_hmm_files(cleaned_ko_list, hmm_directory, hmm_database):
    ko_list_df = pd.read_csv(cleaned_ko_list)
    for knum in ko_list_df['knum']:
        hmm_filepath = os.path.join(hmm_directory, f"{knum}.hmm")
        cmd = f"cat {hmm_filepath} >> {hmm_database}"
        print(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)


def main():
    parser = argparse.ArgumentParser(description="Concatenate multiple HMM profile files into a single HMM database file.")
    parser.add_argument('--cleaned_ko_list', type=str, help="Path to the cleaned KO list.")
    parser.add_argument('--hmm_directory', type=str, help="Directory containing the HMM profile files.")
    parser.add_argument('--hmm_database', type=str, help="Path to the output HMM database file.")
    args = parser.parse_args()

    concatenate_hmm_files(args.cleaned_ko_list, args.hmm_directory, args.hmm_database)


if __name__ == "__main__":
    main()
