#!/bin/bash 
#SBATCH --job-name=37_mcl_preparation
#SBATCH --account=pupko-users
#SBATCH --partition=pupko
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_6_mcl_input_files/37_mcl_preparation_%j.out
#SBATCH --error=/groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_6_mcl_input_files/37_mcl_preparation_%j.err
echo Job ID: $SLURM_JOB_ID
echo Running on nodes: $SLURM_JOB_NODELIST
echo Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE
echo Memory per node: $SLURM_MEM_PER_NODE MB
echo Job name: $SLURM_JOB_NAME
/groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_6_mcl_input_files/37_mcl_preparation.sh
scontrol show job $SLURM_JOB_ID | grep -m 1 "RunTime=" >> /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_6_mcl_input_files/$(echo $SLURM_JOB_NAME)_$(echo $SLURM_JOB_ID)_log.txt
scontrol show job $SLURM_JOB_ID | grep -m 1 "NumCPUs=" >> /groups/pupko/yairshimony/microbializer_runs/5_genomes/outputs/tmp/05_6_mcl_input_files/$(echo $SLURM_JOB_NAME)_$(echo $SLURM_JOB_ID)_log.txt
