import subprocess
import argparse


def main():
    parser = argparse.ArgumentParser(description="Compress and index the concatenated HMM database using hmmpress.")
    parser.add_argument('--hmm_database', type=str, help="Path to the concatenated HMM database file.")

    args = parser.parse_args()

    subprocess.run(["hmmpress", args.hmm_database])


if __name__ == "__main__":
    main()
