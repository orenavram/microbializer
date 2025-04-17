from pathlib import Path
import sys
import logging

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.steps.reconstruct_species_phylogeny import draw_tree
from pipeline.auxiliaries.consts import PROJECT_ROOT_DIR

NEWICK_PATH_UNROOTED = PROJECT_ROOT_DIR / 'gallery' / 'chlamydia_run_b' / 'M1CR0B1AL1Z3R_outputs' / '09_species_phylogeny' / 'final_species_tree.newick'
NEWICK_PATH_ROOTED = PROJECT_ROOT_DIR / 'gallery' / 'chlamydia_run_a' / 'M1CR0B1AL1Z3R_outputs' / '09_species_phylogeny' / 'final_species_tree.newick'


def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    draw_tree(logger, NEWICK_PATH_UNROOTED, False, None)
    draw_tree(logger, NEWICK_PATH_ROOTED, True, "Waddlia_chondrophila_WSU_86-1044")


if __name__ == "__main__":
    main()
