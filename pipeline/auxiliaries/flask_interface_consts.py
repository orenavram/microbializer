
WEBSERVER_NAME = 'M1CR0B1AL1Z3R'

# Arguments keys to run the pipeline with
ARGS_JSON_PATH_KEY = "args_json_path"
CONTIGS_DIR = "contigs_dir"  # can be a path of a directory or a zipped file
IDENTITY_CUTOFF = "identity_cutoff"
E_VALUE_CUTOFF = "e_value_cutoff"
CORE_MINIMAL_PERCENTAGE = "core_minimal_percentage"
BOOTSTRAP = "bootstrap"
OUTGROUP = "outgroup"
FILTER_OUT_PLASMIDS = "filter_out_plasmids"
INPUTS_FASTAS_TYPE = "inputs_fastas_type"
FILE_NAME_FOR_ARGUMENTS = "arguments.json"

# Error description file (the path is relative to the unique folder of the job)
ERROR_FILE_PATH = "error.txt"

# Output files (the paths are relative to the unique folder of the job)
ALL_OUTPUTS_DIRECTORY = WEBSERVER_NAME + "_outputs"
ALL_OUTPUTS_ZIPPED_FORMAT = WEBSERVER_NAME + "_outputs.zip"

ORFS_COUNT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_counts.json"
ORFS_COUNT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_counts.png"
GC_CONTENT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_gc_contents.json"
GC_CONTENT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_gc_contents.png"

OG_TABLE = f"{ALL_OUTPUTS_DIRECTORY}/11_final_table/final_orthologs_table.csv"
PHYLETIC_PATTERN = f"{ALL_OUTPUTS_DIRECTORY}/11_final_table/phyletic_pattern.fas"
OG_SIZE_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/19_groups_sizes_frequency/groups_sizes_frequency.png"

SPECIES_TREE = f"{ALL_OUTPUTS_DIRECTORY}/16_species_phylogeny/final_species_tree.txt"

# Microbializer processor Job variables
MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME = 'lifesciweb'
NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB = '1'
MICROBIALIZER_PROCESSOR_JOB_PREFIX = 'MC'
MICROBIALIZER_PROCESSOR_RESULTS_FILE_NAME = 'results.txt'


MICROBIALIZER_JOB_TEMPLATE = '''
#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q {queue_name}
#PBS -l ncpus={cpu_number}
#PBS -l mem={mem_req}gb
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N {job_name}
#PBS -e {error_files_path}
#PBS -o {output_files_path}

source /powerapps/share/miniconda3-4.7.12/etc/profile.d/conda.sh
conda activate microbializer
pwd
sleep {sleep_interval}

python "/bioseq/microbializer_v2/pipeline/main.py" --{args_json_path_key} {args_json_path}
cat "OKAY" > {results_file_path}

'''
