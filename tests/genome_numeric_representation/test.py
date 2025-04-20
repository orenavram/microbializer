from pathlib import Path
import sys
import logging

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.steps.genome_numeric_representation import get_genome_numeric_representation

GALLERY_RESULTS = SCRIPT_DIR.parent.parent / 'gallery' / '5_genomes' / 'M1CR0B1AL1Z3R_outputs'
ORFS_DIR = GALLERY_RESULTS / '02a_orfs'
ORTHOGROUPS_TABLE = GALLERY_RESULTS / '05a_orthogroups' / 'orthogroups.csv'


def main():
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    get_genome_numeric_representation(logger, ORTHOGROUPS_TABLE, ORFS_DIR, SCRIPT_DIR)


if __name__ == '__main__':
    main()
