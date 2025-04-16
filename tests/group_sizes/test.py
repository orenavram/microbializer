from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.logic_utils import count_and_plot_orthogroups_sizes
from pipeline.auxiliaries.consts import PROJECT_ROOT_DIR

OG_PATH = PROJECT_ROOT_DIR / 'gallery' / 'chlamydia_run_b'  / 'M1CR0B1AL1Z3R_outputs' / '05a_orthogroups' / 'orthogroups.csv'


def main():
    count_and_plot_orthogroups_sizes(OG_PATH, SCRIPT_DIR)


if __name__ == "__main__":
    main()
