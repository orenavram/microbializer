import os
import sys
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from auxiliaries.q_submitter_power import submit_cmds_from_file_to_q

logger = logging.getLogger()

new_line_delimiter = '!@#'

GENOMES_DIR = r"/groups/pupko/yairshimony/test/fastani/AllGenomes"
GENOMES_LIST_PATH = r"/groups/pupko/yairshimony/test/fastani/genomes_paths.txt"
OUTPUTS_DIR = r"/groups/pupko/yairshimony/test/fastani/outputs"
TMP_DIR = r"/groups/pupko/yairshimony/test/fastani/tmp"

# for file_name in os.listdir(GENOMES_DIR):
#     if '.gz' in file_name:
#         file_path = os.path.join(GENOMES_DIR, file_name)
#         cmd = f'gzip -d {file_path}'
#         subprocess.run(cmd, shell=True)



# genomes_paths = [os.path.join(GENOMES_DIR, genome_file_name) for genome_file_name in os.listdir(GENOMES_DIR)]
# with open(GENOMES_LIST_PATH, 'w') as genomes_list_file:
#     genomes_list_file.write('\n'.join(genomes_paths))


for genome_file_name in os.listdir(GENOMES_DIR):
    shell_cmds_as_str = ''
    shell_cmds_as_str += f'source ~/.bashrc{new_line_delimiter}'
    shell_cmds_as_str += f'conda activate ani{new_line_delimiter}'
    shell_cmds_as_str += f'export PATH=$CONDA_PREFIX/bin:$PATH{new_line_delimiter}'
    genome_file_path = os.path.join(GENOMES_DIR, genome_file_name)
    output_path = f'{os.path.join(OUTPUTS_DIR, os.path.splitext(genome_file_name)[0])}_to_all.tsv'
    shell_cmds_as_str += f'fastANI -q {genome_file_path} --rl {GENOMES_LIST_PATH} -o {output_path}{new_line_delimiter}'
    cmds_path = os.path.join(TMP_DIR, f'{genome_file_name}.cmds')
    with open(cmds_path, 'w') as f:
        f.write(f'{shell_cmds_as_str}\t{genome_file_name}\n')  # ADDING THE JOB NAME

    submit_cmds_from_file_to_q(logger, cmds_path, TMP_DIR, 'power-pupko', str(1))

