import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from auxiliaries.pipeline_auxiliaries import get_recursive_step_cummulative_times

PATH = r"C:\Users\yairs\Downloads\temp"


def main():
    get_recursive_step_cummulative_times(PATH)


if __name__ == "__main__":
    main()
