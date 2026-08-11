#!/bin/bash

#SBATCH --job-name=weekly_report
#SBATCH --account=pupko-users_v2
#SBATCH --partition=pupko-pool
#SBATCH --qos=owner
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1:00:00
#SBATCH --output=/microbializer/weekly_report_logs/weekly_report_%j.out
#SBATCH --error=/microbializer/weekly_report_logs/weekly_report_%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/microbializer
source /microbializer/miniconda3/etc/profile.d/conda.sh
conda activate /microbializer/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python /microbializer/scripts/weekly_report.py
