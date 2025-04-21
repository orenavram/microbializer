#!/bin/bash

#SBATCH --job-name=microbializer
#SBATCH --account=pupko=users
#SBATCH --partition=pupko
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192
#SBATCH --time=12:00:00
#SBATCH --output=/groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/main_%j.out
#SBATCH --error=/groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/outputs/main_%j.err


echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"

export HOME=/groups/pupko/yairshimony
source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

echo "PATH: $PATH"

python ~/microbializer/pipeline/main.py --args_json_path /groups/pupko/yairshimony/microbializer_runs/chlamydia_run_a/input_parameters.json
