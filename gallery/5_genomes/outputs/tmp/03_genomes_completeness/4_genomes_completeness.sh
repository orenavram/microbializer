source /home/ai_center/ai_users/yairshimony/miniconda/etc/profile.d/conda.sh
conda activate /home/ai_center/ai_users/yairshimony/miniconda/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /home/ai_center/ai_users/yairshimony/microbializer/pipeline/steps/assessing_genome_completeness.py /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/03_genomes_completeness/job_inputs/4.txt /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/03_genomes_completeness/individual_proteomes_outputs -v False --logs_dir /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/03_genomes_completeness --error_file_path /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/M1CR0B1AL1Z3R_outputs/error.txt --job_name 4_genomes_completeness --use_job_manager True --cpus 1
touch /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/03_genomes_completeness/4_genomes_completeness.done
