source /home/ai_center/ai_users/yairshimony/miniconda/etc/profile.d/conda.sh
conda activate /home/ai_center/ai_users/yairshimony/miniconda/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH
export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/11_2_species_phylogeny/xdg_runtime_dir
python /home/ai_center/ai_users/yairshimony/microbializer/pipeline/steps/reconstruct_species_phylogeny.py /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/steps_results/09_3_aligned_core_proteome_reduced/aligned_core_proteome.fas /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/steps_results/11_2_species_phylogeny/final_species_tree.newick /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/11_2_species_phylogeny --bootstrap True -v False --logs_dir /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/11_2_species_phylogeny --error_file_path /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/M1CR0B1AL1Z3R_outputs/error.txt --job_name tree_reconstruction --use_job_manager True --cpus 20
touch /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/tmp/11_2_species_phylogeny/tree_reconstruction.done
