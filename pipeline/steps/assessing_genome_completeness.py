from sys import argv
import os
import subprocess
import argparse
import logging
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.pipeline_auxiliaries import get_job_logger


BACTERIA_CORE_GENES_HMM_PROFILES_PATH = '/groups/pupko/naamawagner/Microbializer/Busco/hmms'
CORE_GENES_COUNT = len(os.listdir(BACTERIA_CORE_GENES_HMM_PROFILES_PATH))


def compute_genome_completeness(genomic_translated_f, out_dir, logger):
    """
    input:
        genomic_translated_f - protein fasta file of one genome
        out_dir - directory to which the results will be written to
    function specification:
        run hmmserach of all the hmm_profiles vs the genomic protein fasta, and based on the results gives the genome completeness score (in percentage of different profiles found in the genome).
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    score = 0
    for profile in os.listdir(BACTERIA_CORE_GENES_HMM_PROFILES_PATH):
        profile_path = os.path.join(BACTERIA_CORE_GENES_HMM_PROFILES_PATH, profile)
        hmmsearch_out_file_path = os.path.join(out_dir, f'{profile.split(".")[0]}.txt')
        cmd = f'hmmsearch -E 0.1 --pfamtblout {hmmsearch_out_file_path} {profile_path} {genomic_translated_f}'
        subprocess.check_output(cmd, shell=True)

        with open(hmmsearch_out_file_path) as out_hmmsearch:
            # Examine the first sequence hit (=the most significant hit = the first line that doesn't start with #)
            for line in out_hmmsearch:
                if not line.startswith('#'):
                    if float(line.split()[2]) < 10 ** (-4):
                        score += 1
                    else:
                        logger.info(f"Proteome {genomic_translated_f} doesn't include a gene that matches the profile {profile}")
                    break

    return round((score / CORE_GENES_COUNT) * 100)


def main(proteome_path, output_dir, logger):
    """
    the main function that computes the genome completeness of the proteome and saves the results to the output dir.
    """
    genome_name = os.path.basename(proteome_path).split('.')[0]
    genome_out_dir = os.path.join(output_dir, genome_name)
    completeness_score = compute_genome_completeness(proteome_path, genome_out_dir, logger)
    with open(os.path.join(genome_out_dir, 'result.txt'), 'w') as fp:
        fp.write(str(completeness_score))


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('proteome_path', help='path to a protein fasta file')
    parser.add_argument('output_dir', help='path to the output dir')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)

    try:
        main(args.proteome_path, args.output_dir, logger)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')
