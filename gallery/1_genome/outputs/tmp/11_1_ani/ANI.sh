source /home/ai_center/ai_users/yairshimony/miniconda/etc/profile.d/conda.sh
conda activate /home/ai_center/ai_users/yairshimony/miniconda/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /home/ai_center/ai_users/yairshimony/microbializer/pipeline/steps/ani.py /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani/temp_results/genomes_list.txt /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/steps_results/11_1_ani -v False --logs_dir /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani --error_file_path /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/M1CR0B1AL1Z3R_outputs/error.txt --job_name ANI --use_job_manager True --cpus 20
touch /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani/ANI.done
