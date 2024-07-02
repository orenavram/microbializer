#!/bin/bash

#SBATCH --job-name=4_genomes
#SBATCH --account=power-general-users
#SBATCH --partition=power-general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/groups/pupko/yairshimony/microbializer_runs/4genomes/output_slurm/4genomes_%j.out
#SBATCH --error=/groups/pupko/yairshimony/microbializer_runs/4genomes/output_slurm/4genomes_%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python ~/microbializer/pipeline/main.py --contigs_dir ~/microbializer_runs/4genomes/inputs --output_dir output_slurm
