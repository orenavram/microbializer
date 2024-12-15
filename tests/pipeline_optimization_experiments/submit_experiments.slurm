#!/bin/bash

#SBATCH --job-name=micro_expirements
#SBATCH --account=pupko-users
#SBATCH --partition=pupko
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=120:00:00
#SBATCH --output=/groups/pupko/yairshimony/microbializer_runs/73_ecoli/%j.out
#SBATCH --error=/groups/pupko/yairshimony/microbializer_runs/73_ecoli/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/groups/pupko/yairshimony

python ~/microbializer/pipeline/tests/pipeline_optimization_experiments/submit_experiments.py --inputs_dir_path /groups/pupko/yairshimony/microbializer_runs/73_ecoli/inputs --partition_name pupko
