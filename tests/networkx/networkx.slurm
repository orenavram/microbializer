#!/bin/bash

#SBATCH --job-name=networkx
#SBATCH --account=gpu-research
#SBATCH --partition=cpu-killable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --output=/home/ai_center/ai_users/yairshimony/microbializer/tests/networkx/%j.out
#SBATCH --error=/home/ai_center/ai_users/yairshimony/microbializer/tests/networkx/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/home/ai_center/ai_users/yairshimony
source ~/miniconda/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits no --parallelize_write_ogs no
#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits no --parallelize_write_ogs threads
#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits no --parallelize_write_ogs processes
#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits threads --parallelize_write_ogs no
#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits threads --parallelize_write_ogs threads
#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits threads --parallelize_write_ogs processes
#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits processes --parallelize_write_ogs no
#python ~/microbializer/tests/networkx/run_clustering.py --num_workers 20 --parallelize_load_hits processes --parallelize_write_ogs threads
python ~/microbializer/tests/networkx/run_clustering.py --num_workers $SLURM_JOB_CPUS_PER_NODE --parallelize_load_hits processes --parallelize_write_ogs processes

#python ~/microbializer/tests/networkx/combine_results.py