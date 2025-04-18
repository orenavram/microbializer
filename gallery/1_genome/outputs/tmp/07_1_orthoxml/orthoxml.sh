source /home/ai_center/ai_users/yairshimony/miniconda/etc/profile.d/conda.sh
conda activate /home/ai_center/ai_users/yairshimony/miniconda/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /home/ai_center/ai_users/yairshimony/microbializer/pipeline/steps/create_orthoxml.py /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/steps_results/05_10_orthogroups_final/orthogroups.csv /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/steps_results/07_1_orthoxml --qfo_benchmark False -v False --logs_dir /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/07_1_orthoxml --error_file_path /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/M1CR0B1AL1Z3R_outputs/error.txt --job_name orthoxml --use_job_manager True --cpus 1
touch /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/07_1_orthoxml/orthoxml.done
