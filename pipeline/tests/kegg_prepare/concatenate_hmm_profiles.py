import subprocess
import os
import argparse


def concatenate_hmm_files(hmm_list_file, hmm_directory, hmm_database):
    with open(hmm_list_file, 'r') as file:
        for line in file:
            hmm_filename = line.strip()
            hmm_filepath = os.path.join(hmm_directory, hmm_filename)
            subprocess.run(f"cat {hmm_filepath} >> {hmm_database}", shell=True)


def main():
    parser = argparse.ArgumentParser(description="Concatenate multiple HMM profile files into a single HMM database file.")
    parser.add_argument('--hmm_list_file', type=str, help="Path to the text file containing list of HMM profile filenames.")
    parser.add_argument('--hmm_directory', type=str, help="Directory containing the HMM profile files.")
    parser.add_argument('--hmm_database', type=str, help="Path to the output HMM database file.")
    args = parser.parse_args()

    concatenate_hmm_files(args.hmm_list_file, args.hmm_directory, args.hmm_database)


if __name__ == "__main__":
    main()
