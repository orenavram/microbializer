from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.main_utils import define_intervals


def main():
    batch_size = 8

    genomes_batches = define_intervals(18, batch_size)
    print('\n'.join([str(interval) for interval in genomes_batches]))


if __name__ == '__main__':
    main()
