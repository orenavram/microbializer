#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --account=pupko-users
#SBATCH --partition=pupko
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=120:00:00
#SBATCH --output={path_to_output_dir}/%j.out
#SBATCH --error={path_to_output_dir}/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/groups/pupko/yairshimony
source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python ~/microbializer_prod/pipeline/main.py --contigs_dir {inputs_dir/inputs_zip} --account_name pupko-users --queue_name pupko --inputs_fasta_type {genomes/orfs}