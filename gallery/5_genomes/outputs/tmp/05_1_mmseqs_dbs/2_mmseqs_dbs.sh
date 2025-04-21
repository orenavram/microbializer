source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/mmseqs_v4/mmseqs2_create_db.py /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_1_mmseqs_dbs/job_inputs/2.txt /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/steps_results/02_1_orfs/orfs_translated /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/steps_results/05_1_mmseqs_dbs -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_1_mmseqs_dbs --error_file_path /groups/pupko/yairshimony/microbializer_runs/5_genomes/M1CR0B1AL1Z3R_outputs/error.txt --job_name 2_mmseqs_dbs --use_job_manager True --cpus 1
touch /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_1_mmseqs_dbs/2_mmseqs_dbs.done
