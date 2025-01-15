#!/bin/bash

#SBATCH --job-name=networkx
#SBATCH --account=pupkoweb-users
#SBATCH --partition=pupkoweb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --output=/lsweb/pupko/microbializer/tests/networkx/%j.out
#SBATCH --error=/lsweb/pupko/microbializer/tests/networkx/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/lsweb/pupko/microbializer
source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python ~/tests/networkx/run_clustering.py --inputs_hits_dir /lsweb/pupko/microbializer/microbializer_runs/salmonella_300/outputs_basic/steps_results/05_4_normalize_scores/ --num_workers $SLURM_JOB_CPUS_PER_NODE --parallelize_load_hits processes --parallelize_write_ogs processes --verbose
