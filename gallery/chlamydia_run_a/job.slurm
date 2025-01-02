#!/bin/bash

#SBATCH --job-name=microbializer
#SBATCH --account=pupkoweb-users
#SBATCH --partition=pupkoweb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4096
#SBATCH --time=168:00:00
#SBATCH --output=/lsweb/pupko/microbializer/user_results/1c60cfef-e466-4987-9565-e06ddac46b82/main_%j.out
#SBATCH --error=/lsweb/pupko/microbializer/user_results/1c60cfef-e466-4987-9565-e06ddac46b82/main_%j.err


sleep 120

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"

export HOME=/lsweb/pupko/microbializer
source /lsweb/pupko/microbializer/miniconda3/etc/profile.d/conda.sh
conda activate /lsweb/pupko/microbializer/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

echo "PATH: $PATH"

python "/lsweb/pupko/microbializer/pipeline/main.py" --args_json_path /lsweb/pupko/microbializer/user_results/1c60cfef-e466-4987-9565-e06ddac46b82/input_parameters.json --account_name pupkoweb-users --queue_name pupkoweb
