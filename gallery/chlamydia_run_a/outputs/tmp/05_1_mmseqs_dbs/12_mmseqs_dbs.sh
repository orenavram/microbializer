source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/mmseqs_v4/mmseqs2_create_db.py /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_1_mmseqs_dbs/job_inputs/12.txt /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/steps_results/02_1_orfs/orfs_translated /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/steps_results/05_1_mmseqs_dbs -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_1_mmseqs_dbs --error_file_path /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/M1CR0B1AL1Z3R_outputs/error.txt --job_name 12_mmseqs_dbs --use_job_manager True --cpus 1
touch /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/05_1_mmseqs_dbs/12_mmseqs_dbs.done
