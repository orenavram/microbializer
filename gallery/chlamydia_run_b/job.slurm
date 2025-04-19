#!/bin/bash

#SBATCH --job-name=microbializer
#SBATCH --account=gpu-research
#SBATCH --partition=killable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192
#SBATCH --time=12:00:00
#SBATCH --output=/home/ai_center/ai_users/yairshimony/microbializer/gallery/chlamydia_run_b/main_%j.out
#SBATCH --error=/home/ai_center/ai_users/yairshimony/microbializer/gallery/chlamydia_run_b/main_%j.err


echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"

export HOME=/home/ai_center/ai_users/yairshimony
source ~/miniconda/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

echo "PATH: $PATH"

python ~/microbializer/pipeline/main.py --args_json_path ~/microbializer/gallery/chlamydia_run_b/input_parameters.json
