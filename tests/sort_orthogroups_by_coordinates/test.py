from pathlib import Path
import sys
import logging

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.logic_utils import sort_orthogroups_df_and_rename_ogs

ORFS_COORDINATES_DIR = SCRIPT_DIR / 'orfs_coordinates'
ORTHOGROUPS_TABLE = SCRIPT_DIR / 'orthogroups.csv'


def main():
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    sort_orthogroups_df_and_rename_ogs(logger, ORTHOGROUPS_TABLE, ORFS_COORDINATES_DIR, SCRIPT_DIR / 'sorted_orthogroups.csv')


if __name__ == '__main__':
    main()
