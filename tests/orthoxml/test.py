import os
import pandas as pd
import logging

from pipeline.steps.orthologs_table_variations import build_orthoxml_and_tsv_output


ORTHOGROUPS_PATH = r"C:\Users\TalPNB22\Downloads\orthogroups.csv"
SCRIPT_DIR = os.path.dirname(__file__)


def main():
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)
    df = pd.read_csv(ORTHOGROUPS_PATH)
    build_orthoxml_and_tsv_output(logger, df, SCRIPT_DIR)


if __name__ == '__main__':
    main()
