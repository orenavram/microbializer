#!/bin/bash

#SBATCH --job-name=test_memory
#SBATCH --account=pupkoweb-users
#SBATCH --partition=pupkoweb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=120:00:00
#SBATCH --output=/lsweb/pupko/microbializer/tests/extract_orfs_memory_test/%j.out
#SBATCH --error=/lsweb/pupko/microbializer/tests/extract_orfs_memory_test/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/groups/pupko/yairshimony
source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python /lsweb/pupko/microbializer/tests/extract_orfs_memory_test/test_memory_usage.py "/lsweb/pupko/microbializer/microbializer_runs/1000_genomes/outputs_approximate/steps_results/02_3_concat_orfs/all_orfs.fna" "/lsweb/pupko/microbializer/microbializer_runs/1000_genomes/outputs_approximate/steps_results/02_3_concat_orfs/all_proteomes.faa"