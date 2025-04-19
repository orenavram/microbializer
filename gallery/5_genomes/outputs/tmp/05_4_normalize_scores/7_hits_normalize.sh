source /home/ai_center/ai_users/yairshimony/miniconda/etc/profile.d/conda.sh
conda activate /home/ai_center/ai_users/yairshimony/miniconda/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /home/ai_center/ai_users/yairshimony/microbializer/pipeline/steps/normalize_hits_scores.py /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_4_normalize_scores/jobs_inputs/7.txt /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/steps_results/05_4_normalize_scores --use_parquet True -v False --logs_dir /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_4_normalize_scores --error_file_path /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/M1CR0B1AL1Z3R_outputs/error.txt --job_name 7_hits_normalize --use_job_manager True --cpus 1
touch /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_4_normalize_scores/7_hits_normalize.done
