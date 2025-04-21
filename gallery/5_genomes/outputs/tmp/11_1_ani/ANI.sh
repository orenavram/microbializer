source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/ani.py /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/11_1_ani/temp_results/genomes_list.txt /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/steps_results/11_1_ani -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/11_1_ani --error_file_path /groups/pupko/yairshimony/microbializer_runs/5_genomes/M1CR0B1AL1Z3R_outputs/error.txt --job_name ANI --use_job_manager True --cpus 20
touch /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/11_1_ani/ANI.done
