#!/bin/bash

#SBATCH --job-name=microbializer
#SBATCH --account=gpu-research
#SBATCH --partition=cpu-killable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=120:00:00
#SBATCH --output=/home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/%j.out
#SBATCH --error=/home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/outputs/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/home/ai_center/ai_users/yairshimony
source ~/miniconda/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python ~/microbializer/pipeline/main.py --args_json_path /home/ai_center/ai_users/yairshimony/microbializer_runs/5_genomes/input_parameters.json