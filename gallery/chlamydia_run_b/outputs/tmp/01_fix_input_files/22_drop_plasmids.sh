source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
python /groups/pupko/yairshimony/microbializer/pipeline/steps/drop_plasmids_and_fix_frames.py /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_b/outputs/tmp/01_fix_input_files/job_inputs/22.txt /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_b/outputs/steps_results/01_fix_input_files --drop_plasmids True --fix_frames False -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_b/outputs/tmp/01_fix_input_files --error_file_path /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_b/M1CR0B1AL1Z3R_outputs/error.txt --job_name 22_drop_plasmids --use_job_manager True --cpus 1
touch /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_b/outputs/tmp/01_fix_input_files/22_drop_plasmids.done
