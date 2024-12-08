import argparse
import os
from pathlib import Path
import subprocess

JOB_TEMPLATE = """
#!/bin/bash

#SBATCH --job-name=microbializer
#SBATCH --account={partition_name}-users
#SBATCH --partition={partition_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={memory}
#SBATCH --output={output_dir_path}/%j.out
#SBATCH --error={output_dir_path}/%j.err

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Job name: $SLURM_JOB_NAME"

export HOME=/groups/pupko/yairshimony
source ~/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

python {microbializer_dir}/pipeline/main.py --contigs_dir {inputs_dir_path} --output_dir {outputs_dir_name} --inputs_fasta_type genomes --account_name {partition_name}-users -q {partition_name} --add_orphan_genes_to_ogs --only_calc_ogs {pre_cluster_flags} {run_optimized_mmseqs} {use_parquet}
"""


MICROBIALIZER_DIR = '/groups/pupko/yairshimony/microbializer'
MAIN_JOB_MEMORY = '4G'

PRE_CLUSTER_FLAGS = {'no': '', 'mergeAtEnd': '--pre_cluster_orthogroups_inference', 'mergeAfterMmseqs': '--pre_cluster_orthogroups_inference --unify_clusters_after_mmseqs'}
OPTIMIZED_MMSEQS_FLAGS = {'no': '', 'yes': '--run_optimized_mmseqs'}
USE_PARQUET_FLAGS = {'no': '', 'yes': '--use_parquet'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--partition_name', help='The partition name', type=str, required=True)
    parser.add_argument('--inputs_dir_path', help='The inputs directory path', type=str, required=True)
    args = parser.parse_args()

    run_dir = Path(args.inputs_dir_path).parent.absolute()

    for pre_cluster_flag_id, pre_cluster_flag in PRE_CLUSTER_FLAGS.items():
        for optimized_mmseqs_flag_id, optimized_mmseqs_flag in OPTIMIZED_MMSEQS_FLAGS.items():
            for use_parquet_flag_id, use_parquet_flag in USE_PARQUET_FLAGS.items():
                flags_str = f'preCluster_{pre_cluster_flag_id}_optimizedMmseqs_{optimized_mmseqs_flag_id}_useParquet_{use_parquet_flag_id}'
                outputs_dir_name = f'outputs_{flags_str}'
                job_content = JOB_TEMPLATE.format(partition_name=args.partition_name,
                                                  memory=MAIN_JOB_MEMORY, output_dir_path=os.path.join(run_dir, outputs_dir_name),
                                                  microbializer_dir=MICROBIALIZER_DIR, inputs_dir_path=args.inputs_dir_path,
                                                  outputs_dir_name=outputs_dir_name, pre_cluster_flags=pre_cluster_flag,
                                                  run_optimized_mmseqs=optimized_mmseqs_flag, use_parquet=use_parquet_flag)

                job_path = f'{run_dir}/job_{flags_str}.slurm'
                with open(job_path, 'w') as job_fp:
                    job_fp.write(job_content)

                subprocess.run(f'sbatch {job_path}', shell=True, capture_output=True, text=True, check=True)


if __name__ == '__main__':
    main()
