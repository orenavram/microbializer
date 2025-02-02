import os
import Bio.Seq
from sys import argv
import argparse
import sys
import logging

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import get_job_logger


def wait_for_output_folder(logger, output_folder, max_waiting_time=300):
    i = 0
    while not os.path.exists(output_folder):
        logger.info(f'Waiting to {output_folder} to be generated... (waited {i} seconds)')
        i += 1
        if i > max_waiting_time:
            raise OSError(
                f'{output_folder} was not generated after {max_waiting_time} second. Failed to continue the analysis.')
        sleep(1)


def get_genome_sequence(genome_path):
    with open(genome_path) as f:
        f.readline()  # get rid of the header
        genome = f.read().replace('\n', '')

    assert '>' not in genome, f'Extracting ORFs and promoter for {genome_path} was failed. This script handles ' \
                              f'only SINGLE contigs genomes. What is the meaning of "global relative location" when you have a set of ' \
                              f'contigs, right?\n'

    return genome


def get_sequence(genome, genome_len, start, stop, reverse_complement, promoters_length):
    if reverse_complement == '-1':
        stop = stop + promoters_length
        if start < stop:
            sequence = genome[start:min(stop, genome_len)]
        else:
            sequence = genome[start:] + genome[:min(stop, genome_len)]
        if stop > genome_len:  # add bases from the beginning of the genome (closed circle)
            sequence += genome[:stop - genome_len]
        sequence = Bio.Seq.reverse_complement(sequence)

    else:  # not a reverse complement
        start = start - promoters_length
        if start < stop:
            sequence = genome[max(start, 0):stop]
        else:
            sequence = genome[max(start, 0):] + genome[:stop]
        if start < 0:  # add bases from the end of the genome
            sequence = genome[start:] + sequence
    return sequence


def extract_promoters_and_orfs(logger, prodigal_orfs_path, genome_path, promoters_length, output_path):
    genome = get_genome_sequence(genome_path)
    genome_len = len(genome)
    logger.info(f'Genome length is: {genome_len}')

    with open(prodigal_orfs_path) as f:
        result = ''
        for line in f:
            if not line.startswith('>'):
                continue
            logger.debug(line)
            header, start, stop, reverse_complement = line.split(' # ')[:4]
            start, stop = int(start) - 1, int(stop)  # Prodigal start is shifted by one
            sequence = get_sequence(genome, genome_len, start, stop, reverse_complement, promoters_length)
            logger.debug(sequence)
            start += 1
            if reverse_complement == '-1':
                stop += promoters_length
                if stop % genome_len < start:
                    stop = stop % genome_len
            else:
                start -= promoters_length
                if start < 1:
                    start -= 1
            result += f'{header} # {start % (genome_len + 1)} # {stop % (genome_len + 1)} # {reverse_complement}\n{sequence}\n'

    with open(output_path, 'w') as f:
        f.write(result)

    logger.info(f'Promoters and ORFs of {genome_path} were extracted successfully to {output_path}')


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_path',
                        help='A path to a genome from which promoters and orfs should be extracted',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('prodigal_orfs_path',
                        help='A path to a Prodigal output file from which coordinates should be extracted',
                        type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('output_path', help='A path to which the promoters and orfs will be written')
    # type=lambda path: path if os.path.exists(os.path.split(path)[0])
    # else parser.error(f'output folder {os.path.split(path)[0]} does not exist!'))
    parser.add_argument('--promoters_length',
                        help='How much bases upstream to the ORF should be extracted',
                        type=int, default=300)
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    # wait_for_output_folder(os.path.split(args.output_path)[0])

    logger.info(script_run_message)
    try:
        extract_promoters_and_orfs(logger, args.prodigal_orfs_path, args.genome_path, args.promoters_length, args.output_path)
    except Exception as e:
        logger.exception(f'Error in {Path(__file__).name}')
