#!/bin/bash

#SBATCH --job-name=microbi
#SBATCH --account=pupkoweb-users
#SBATCH --partition=pupkoweb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=120:00:00
#SBATCH --output=/lsweb/pupko/microbializer/microbializer_runs/4_genomes/outputs/%j.out
#SBATCH --error=/lsweb/pupko/microbializer/microbializer_runs/4_genomes/outputs/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/lsweb/pupko/microbializer
source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python /groups/pupko/yairshimony/microbializer/pipeline/main.py --contigs_dir /lsweb/pupko/microbializer/microbializer_runs/4_genomes/4_genomes.zip --account_name pupkoweb-users --queue_name pupkoweb --inputs_fasta_type genomes --mmseqs_use_dbs --use_parquet