source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/create_orthoxml.py /groups/pupko/yairshimony/microbializer_runs/1_genome_small/outputs/steps_results/05_11_sort_orthogroups_by_coordinates/orthogroups.csv /groups/pupko/yairshimony/microbializer_runs/1_genome_small/outputs/steps_results/07_1_orthoxml --qfo_benchmark False -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/1_genome_small/outputs/tmp/07_1_orthoxml --error_file_path /groups/pupko/yairshimony/microbializer_runs/1_genome_small/M1CR0B1AL1Z3R_outputs/error.txt --job_name orthoxml --use_job_manager True --cpus 1
touch /groups/pupko/yairshimony/microbializer_runs/1_genome_small/outputs/tmp/07_1_orthoxml/orthoxml.done
