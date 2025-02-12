import argparse
from pathlib import Path
import subprocess
import sys
import shutil

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent.parent))

from pipeline.auxiliaries.run_step_utils import add_default_step_args, run_step


def mcl(logger, input_file, output_file, cpus):
    if output_file.exists():
        return

    # --abc for a columns format, i.e., item1\item2\tscore
    cmd = f'mcl "{input_file}" -I 1.5 --abc -o "{output_file}" -te {cpus} -V all'
    logger.info(f'Starting MCL. Calling: {cmd}')
    subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    logger.info(f'MCL finished. Output written to {output_file}')


def verify(logger, og_name, input_file, output_dir):
    with open(input_file, 'r') as f:
        lines = [line.rstrip() for line in f if line]

    if len(lines) == 0:
        raise ValueError(f'{input_file} is empty! There\'s a bug in the previous step!')
    elif len(lines) == 1:
        output_file_path = output_dir / f'{og_name}.verified_cluster'
        if output_file_path.exists():
            return
        shutil.copyfile(input_file, output_file_path)
        logger.info(f'Only one cluster in {input_file}. Copied it to {output_file_path}')
    else:  # 1 < len(lines)
        og_subset_id = 0
        for line in lines:
            if len(line.split('\t')) == 1:
                continue  # Skip lines with 1 protein (orphan genes)
            verified_cluster_path = output_dir / f"{og_name}_{og_subset_id}.split_cluster"
            if not verified_cluster_path.exists():
                with open(verified_cluster_path, 'w') as verified_cluster_file:
                    verified_cluster_file.write(line)
            og_subset_id += 1
        logger.info(f'{input_file} was split into {og_subset_id} clusters in {output_dir}')


def run_mcl_on_all_putative_ogs(logger, mcl_input_dir, job_input_path, mcl_output_dir, verified_clusters_dir, cpus):
    with open(job_input_path, 'r') as f:
        ogs_names = [line.strip() for line in f]

    for og_name in ogs_names:
        input_mcl_path = mcl_input_dir / f'{og_name}.mcl_input'
        output_mcl_path = mcl_output_dir / f'{og_name}.mcl_output'
        mcl(logger, input_mcl_path, output_mcl_path, cpus)
        verify(logger, og_name, output_mcl_path, verified_clusters_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mcl_input_dir', type=Path, help='path to dir of mcl input files')
    parser.add_argument('job_input_path', type=Path, help='')
    parser.add_argument('mcl_output_dir', type=Path, help='path to dir the MCL analysis will be written')
    parser.add_argument('verified_clusters_dir', type=Path, help='dir path to which verified clusters are written')
    add_default_step_args(parser)
    args = parser.parse_args()

    run_step(args, run_mcl_on_all_putative_ogs, args.mcl_input_dir, args.job_input_path, args.mcl_output_dir,
             args.verified_clusters_dir, args.cpus)
