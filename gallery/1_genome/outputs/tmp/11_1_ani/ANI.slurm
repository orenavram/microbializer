#!/bin/bash 
#SBATCH --job-name=ANI
#SBATCH --account=gpu-research
#SBATCH --partition=cpu-killable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=/home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani/ANI_%j.out
#SBATCH --error=/home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani/ANI_%j.err
echo Job ID: $SLURM_JOB_ID
echo Running on nodes: $SLURM_JOB_NODELIST
echo Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE
echo Memory per node: $SLURM_MEM_PER_NODE MB
echo Job name: $SLURM_JOB_NAME
/home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani/ANI.sh
scontrol show job $SLURM_JOB_ID | grep -m 1 "RunTime=" >> /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani/$(echo $SLURM_JOB_NAME)_$(echo $SLURM_JOB_ID)_log.txt
scontrol show job $SLURM_JOB_ID | grep -m 1 "NumCPUs=" >> /home/ai_center/ai_users/yairshimony/microbializer_runs/1_genome_no_paralogs/outputs/tmp/11_1_ani/$(echo $SLURM_JOB_NAME)_$(echo $SLURM_JOB_ID)_log.txt
