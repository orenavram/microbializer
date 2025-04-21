source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/assessing_genome_completeness.py /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/03_genomes_completeness/job_inputs/10.txt /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/03_genomes_completeness/individual_proteomes_outputs -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/03_genomes_completeness --error_file_path /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/M1CR0B1AL1Z3R_outputs/error.txt --job_name 10_genomes_completeness --use_job_manager True --cpus 1
touch /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/03_genomes_completeness/10_genomes_completeness.done
