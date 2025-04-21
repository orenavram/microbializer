source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/construct_putative_orthologs_table.py /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/steps_results/05_4_normalize_scores /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/steps_results/05_5_putative_table/putative_orthologs_table.csv -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_5_putative_table --error_file_path /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/M1CR0B1AL1Z3R_outputs/error.txt --job_name putative_table --use_job_manager True --cpus 1
touch /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_5_putative_table/putative_table.done
