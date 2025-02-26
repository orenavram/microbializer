import sys
import argparse
from pathlib import Path
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step
from pipeline.auxiliaries.general_utils import str_to_bool


def filter_out_plasmids(logger, input_genome_path, output_genome_path, drop_plasmids, fix_frames):
    """
        input_genome_path: path to an input fasta file with prokaryotic genome
        output_genome_path: path to a filtered genome without plasmids
    """
    if drop_plasmids:
        logger.info(f'Removing plasmids from {input_genome_path}...')
    if fix_frames:
        logger.info(f'Fixing frames of {input_genome_path}...')

    new_records = []
    for record in SeqIO.parse(input_genome_path, 'fasta'):
        if drop_plasmids and ('plasmid' in record.id.lower() or 'plasmid' in record.description.lower()):
            logger.info(f'Dropped plasmid sequence {record.id}')
        else:
            new_records.append(record)

        if fix_frames:
            if 'frame=2' in record.description:
                record.seq = record.seq[1:]
                logger.info(f'Fixed frame of {record.id} (frame=2)')
            elif 'frame=3' in record.description:
                record.seq = record.seq[2:]
                logger.info(f'Fixed frame of {record.id} (frame=3)')

    if not new_records:
        logger.info(f'No records left for {input_genome_path} (probably contained only plasmids)')
        return

    SeqIO.write(new_records, output_genome_path, 'fasta')


def filter_out_plasmids_of_all_files(logger, job_input_path, output_dir, drop_plasmids, fix_frames):
    with open(job_input_path, 'r') as f:
        for line in f:
            genome_path = Path(line.strip())
            genome_file_name = genome_path.name
            output_path = output_dir / genome_file_name
            filter_out_plasmids(logger, genome_path, output_path, drop_plasmids, fix_frames)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path,
                        help='path to a file that contains the genome names to drop plasmids from')
    parser.add_argument('output_dir', type=Path, help='path to output dir')
    parser.add_argument('--drop_plasmids', type=str_to_bool, help='Drop plasmids from the genome file')
    parser.add_argument('--fix_frames', type=str_to_bool, help='Fix frames of the genome file')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, filter_out_plasmids_of_all_files, args.job_input_path, args.output_dir, args.drop_plasmids,
             args.fix_frames)
