from pathlib import Path
import sys
import pandas as pd
import logging

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.steps.orthologs_table_variations import plot_presence_absence_matrix
from pipeline.auxiliaries.consts import PROJECT_ROOT_DIR

OG_PATH_40 = PROJECT_ROOT_DIR / 'gallery' / 'chlamydia_run_b' / 'M1CR0B1AL1Z3R_outputs' / '05a_orthogroups' / 'orthogroups.csv'
OG_PATH_300 = SCRIPT_DIR / 'final_orthologs_table_300.csv'
OG_PATH_4 = PROJECT_ROOT_DIR / 'gallery' / '4_genomes' / 'M1CR0B1AL1Z3R_outputs' / '05a_orthogroups' / 'orthogroups.csv'


def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    orthogroups_df = pd.read_csv(OG_PATH_4, dtype=str)
    plot_presence_absence_matrix(logger, orthogroups_df, SCRIPT_DIR)


if __name__ == "__main__":
    main()
