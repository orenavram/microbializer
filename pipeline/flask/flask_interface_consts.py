
WEBSERVER_NAME = 'M1CR0B1AL1Z3R'

# Arguments keys to run the pipeline with
ARGS_JSON_PATH_KEY = "args_json_path"
JOB_PARAMETERS_FILE_NAME = "input_parameters.json"
CONTIGS_DIR = "contigs_dir"  # can be a path of a directory or a zipped file

# Input parametrs from users
IDENTITY_CUTOFF = "identity_cutoff"
E_VALUE_CUTOFF = "e_value_cutoff"
CORE_MINIMAL_PERCENTAGE = "core_minimal_percentage"
COVERAGE_CUTOFF = "coverage_cutoff"
BOOTSTRAP = "bootstrap"
OUTGROUP = "outgroup"
FILTER_OUT_PLASMIDS = "filter_out_plasmids"
ADD_ORPHAN_GENES_TO_OGS = "add_orphan_genes_to_ogs"
INPUT_FASTA_TYPE = "inputs_fasta_type"

INPUTS_ARE_ANNOTATED_PROTEOMES = "inputs_are_annotated_proteomes"


FILE_NAME_FOR_ARGUMENTS = "arguments.json"

# Error description file (the path is relative to the unique folder of the job)
ERROR_FILE_PATH = "error.txt"

# Output files (the paths are relative to the unique folder of the job)
ALL_OUTPUTS_DIRECTORY = WEBSERVER_NAME + "_outputs"
ALL_OUTPUTS_ZIPPED_FORMAT = WEBSERVER_NAME + "_outputs.zip"

ORFS_COUNT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plots/orfs_counts.json"
ORFS_COUNT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plots/orfs_counts.png"
GC_CONTENT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plots/orfs_gc_contents.json"
GC_CONTENT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plots/orfs_gc_contents.png"

OG_TABLE = f"{ALL_OUTPUTS_DIRECTORY}/11_final_table/final_orthologs_table.csv"
PHYLETIC_PATTERN = f"{ALL_OUTPUTS_DIRECTORY}/11_final_table/phyletic_pattern.fas"
OG_SIZE_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/19_groups_sizes_frequency/groups_sizes_frequency.png"

SPECIES_TREE = f"{ALL_OUTPUTS_DIRECTORY}/16_species_phylogeny/final_species_tree.txt"

TITLE_HISTORGRAM_FOR_ORFS = "Open Reading Frames (ORFs)"
TITLE_HISTORGRAM_FOR_GC_CONTENT = "GC content"

DATA_2_VIEW_IN_HISTOGRAM = {
    # use title as key and path to json as value
    # REMEMBER! json file should be a dict with the genome name as key and a scalar as the value
    TITLE_HISTORGRAM_FOR_ORFS: ORFS_COUNT_PER_GENOME,
    TITLE_HISTORGRAM_FOR_GC_CONTENT: GC_CONTENT_PER_GENOME
}

# Microbializer processor Job variables
MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME = 'pupkoweb'
NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB = '1'
MICROBIALIZER_PROCESSOR_JOB_PREFIX = 'MC'
MICROBIALIZER_PROCESSOR_RESULTS_FILE_NAME = 'results.txt'


MICROBIALIZER_JOB_TEMPLATE = '''#!/bin/bash

sleep {sleep_interval}

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"

#python "/lsweb/pupko/microbializer/pipeline/main.py" --{args_json_path_key} {args_json_path}
echo OKAY > {results_file_path}

'''
