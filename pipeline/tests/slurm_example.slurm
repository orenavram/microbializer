#!/bin/bash

#SBATCH --job-name=10_genomes
#SBATCH --account=pupko-users
#SBATCH --partition=pupko
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/groups/pupko/yairshimony/microbializer_runs/example_data/outputs/4genomes_%j.out
#SBATCH --error=/groups/pupko/yairshimony/microbializer_runs/example_data/outputs/4genomes_%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"

export HOME=/groups/pupko/yairshimony
source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python ~/microbializer_prod/pipeline/main.py --contigs_dir ~/microbializer_runs/example_data/inputs
