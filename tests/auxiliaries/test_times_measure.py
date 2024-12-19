import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from pipeline.auxiliaries.pipeline_auxiliaries import get_jobs_cummulative_time

PATH = SCRIPT_DIR


def main():
    get_jobs_cummulative_time(PATH)


if __name__ == "__main__":
    main()
