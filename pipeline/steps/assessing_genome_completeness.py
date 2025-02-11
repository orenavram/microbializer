import subprocess
import argparse
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from auxiliaries.pipeline_auxiliaries import add_default_step_args, run_step
from auxiliaries import consts


def compute_genome_completeness(logger, genomic_translated_f, out_dir):
    """
    input:
        genomic_translated_f - protein fasta file of one genome
        out_dir - directory to which the results will be written to
    function specification:
        run hmmserach of all the hmm_profiles vs the genomic protein fasta, and based on the results gives the genome completeness score (in percentage of different profiles found in the genome).
    """
    score = 0
    for profile_path in consts.BACTERIA_CORE_GENES_HMM_PROFILES_PATH.iterdir():
        hmmsearch_out_file_path = out_dir / f'{profile_path.stem}.txt'
        cmd = f'hmmsearch --noali -o /dev/null -E 0.1 --pfamtblout {hmmsearch_out_file_path} {profile_path} {genomic_translated_f}'
        logger.info(f'Running command: {cmd}')
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

        with open(hmmsearch_out_file_path) as out_hmmsearch:
            # Examine the first sequence hit (=the most significant hit = the first line that doesn't start with #)
            for line in out_hmmsearch:
                if not line.startswith('#'):
                    if line.isspace() or float(line.split()[2]) > consts.BUSCO_EVAULE_CUTOFF:
                        logger.info(
                            f"Proteome {genomic_translated_f} doesn't include a gene that matches the profile {profile_path.name}")
                    else:
                        score += 1
                    break

    return (score / consts.BACTERIA_CORE_GENES_COUNT) * 100


def main(logger, job_input_path, output_dir):
    """
    the main function that computes the genome completeness of the proteome and saves the results to the output dir.
    """
    with open(job_input_path, 'r') as f:
        for line in f:
            proteome_path = line.strip()
            strain_name = Path(proteome_path).stem
            strain_out_dir = output_dir / strain_name
            strain_out_dir.mkdir(exist_ok=True)
            completeness_score = compute_genome_completeness(logger, proteome_path, strain_out_dir)
            proteome_score_path = strain_out_dir / 'result.txt'
            with open(proteome_score_path, 'w') as fp:
                fp.write(str(completeness_score))
            logger.info(f'Genome completeness score for {proteome_path} was written to {proteome_score_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('job_input_path', type=Path,
                        help='path to a file that contains the genome names to asses completeness for')
    parser.add_argument('output_dir', type=Path, help='path to the output dir')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, main, args.job_input_path, args.output_dir)
