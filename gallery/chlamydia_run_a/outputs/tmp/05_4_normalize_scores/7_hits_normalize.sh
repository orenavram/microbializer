source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/normalize_hits_scores.py /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_4_normalize_scores/jobs_inputs/7.txt /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/steps_results/05_4_normalize_scores --use_parquet True -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_4_normalize_scores --error_file_path /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/M1CR0B1AL1Z3R_outputs/error.txt --job_name 7_hits_normalize --use_job_manager True --cpus 1
touch /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_4_normalize_scores/7_hits_normalize.done
