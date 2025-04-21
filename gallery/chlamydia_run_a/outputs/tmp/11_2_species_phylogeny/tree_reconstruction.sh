source /groups/pupko/yairshimony/miniconda3/etc/profile.d/conda.sh
conda activate /groups/pupko/yairshimony/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/11_2_species_phylogeny/xdg_runtime_dir
python /groups/pupko/yairshimony/microbializer/pipeline/steps/reconstruct_species_phylogeny.py /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/steps_results/09_3_aligned_core_proteome_reduced/aligned_core_proteome.fas /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/steps_results/11_2_species_phylogeny/final_species_tree.newick /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/11_2_species_phylogeny --bootstrap True --outgroup Waddlia_chondrophila_WSU_86-1044 -v False --logs_dir /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/11_2_species_phylogeny --error_file_path /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/M1CR0B1AL1Z3R_outputs/error.txt --job_name tree_reconstruction --use_job_manager True --cpus 20
touch /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/tmp/11_2_species_phylogeny/tree_reconstruction.done
