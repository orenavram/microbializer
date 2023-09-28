from sys import argv
import os
import subprocess
import argparse
import logging
import csv
import shutil

from auxiliaries.pipeline_auxiliaries import get_job_logger

def compute_genome_completeness(genomic_translated_f,hmm_dir,out_dir):
    '''
    input:
        genomic_translated_f - protein fasta file of one genome
        hmm_dir - directory that contains all the hmm profiles to search (and only hmm profiles)
        out_dir - directory to which the results will be written to
    function specification:
        run hmmserach of all the hmm_profiles vs the genomic protein fasta, and based on the results gives the genome completeness score (in percentage of different profiles found in the genome).
    '''
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    score = 0
    desired_score = 0
    for profile in os.listdir(hmm_dir):
        cmd = f'hmmsearch -E 1 --pfamtblout {out_dir}/{profile.split(".")[0]}.txt {hmm_dir}/{profile} {genomic_translated_f}'
        subprocess.check_output(cmd,shell=True)

        desired_score += 1
        with open(os.path.join(out_dir,f'{profile.split(".")[0]}.txt')) as out_hmmsearch:
            for line in out_hmmsearch:
                if not line.startswith('#'):
                    if float(line.split()[2]) < 10**(-4):
                        score += 1
                        break
    return round((score/desired_score)*100)
        
def main(translated_genomes_dir,wd,hmm_profiles_dir='/groups/pupko/naamawagner/Microbializer/Busco/hmms'):
    '''
    the main function that computes the genome completeness of all the genomes and saves the results to the output file.
    input specification:
        translated_genomes_dir - directory in which all the protein fasta files that represent the different genomes are found
        hmm_out_dir - a directory to which the hmmerseacrh results will be created
        out_f - a csv file that will contain the scores of the different genomes'
        hmm_profiles_dir = directory that contains all the hmm profiles to search (and only hmm profiles)
    '''
    os.chdir(wd)
    hmm_out_dir = os.path.join(wd,'out_busco_hmm_search')
    out_f = 'genomes_completeness_assessment.csv'
    completeness_dict = {}
    if not os.path.exists(hmm_out_dir):
        os.makedirs(hmm_out_dir)
    for prot_f in os.listdir(translated_genomes_dir):
        genome_name = prot_f.split('.')[0]
        prot_f_path = os.path.join(translated_genomes_dir,prot_f)
        out_dir = os.path.join(hmm_out_dir,genome_name)
        complete_score = compute_genome_completeness(prot_f_path,hmm_profiles_dir,out_dir)
        completeness_dict[genome_name]=complete_score
    with open(out_f,'w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['genome','completeness_percentage'])
        for genome in completeness_dict:
            writer.writerow([genome,completeness_dict[genome]])
    # comment the next line if you don't wish to delete hmmer results
    shutil.rmtree(hmm_out_dir)


if __name__ == '__main__':
    script_run_message = f'Starting command is: {" ".join(argv)}'
    #print(script_run_message)

    parser = argparse.ArgumentParser()
    parser.add_argument('translated_genomes_dir', help='path to a directory with protein fasta files of the different genomes')
    parser.add_argument('working_dir', help='path to the working directory')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    parser.add_argument('--logs_dir', help='path to tmp dir to write logs to')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_job_logger(args.logs_dir, level)

    logger.info(script_run_message)

    try:
        main(args.translated_genomes_dir, args.working_dir)
    except Exception as e:
        logger.exception(f'Error in {os.path.basename(__file__)}')