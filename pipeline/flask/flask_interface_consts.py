
WEBSERVER_NAME = 'M1CR0B1AL1Z3R'

# Arguments keys to run the pipeline with
ARGS_JSON_PATH_KEY = "args_json_path"
JOB_PARAMETERS_FILE_NAME = "input_parameters.json"
RUN_DIR = "run_dir"

# Input parameters from users
IDENTITY_CUTOFF = "identity_cutoff"
E_VALUE_CUTOFF = "e_value_cutoff"
CORE_MINIMAL_PERCENTAGE = "core_minimal_percentage"
COVERAGE_CUTOFF = "coverage_cutoff"
BOOTSTRAP = "bootstrap"
OUTGROUP = "outgroup"
FILTER_OUT_PLASMIDS = "filter_out_plasmids"
ADD_ORPHAN_GENES_TO_OGS = "add_orphan_genes_to_ogs"
INPUT_FASTA_TYPE = "inputs_fasta_type"

# Error description file (the path is relative to the unique folder of the job)
ERROR_FILE_PATH = "error.txt"

# Output files (the paths are relative to the unique folder of the job)
ALL_OUTPUTS_DIRECTORY = WEBSERVER_NAME + "_outputs"
ALL_OUTPUTS_ZIPPED_FORMAT = WEBSERVER_NAME + "_outputs.zip"

ANI_CSV = f"{ALL_OUTPUTS_DIRECTORY}/01_ANI/ani_pairwise_values.csv"
ANI_MAP = f"{ALL_OUTPUTS_DIRECTORY}/01_ANI/ani_map.json"

ORFS_COUNT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/02b_orfs_plots/orfs_counts.json"
ORFS_COUNT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/02b_orfs_plots/orfs_counts.png"
GC_CONTENT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/02b_orfs_plots/orfs_gc_content.json"
GC_CONTENT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/02b_orfs_plots/orfs_gc_content.png"
GENOME_COMPLETENESS_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/03_genomes_completeness/genomes_completeness.json"
GENOME_COMPLETENESS_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/03_genomes_completeness/genomes_completeness.png"
ORPHAN_GENES_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/04_orphan_genes/orphan_genes_count.json"
ORPHAN_GENES_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/04_orphan_genes/orphan_genes_count.png"

OG_TABLE = f"{ALL_OUTPUTS_DIRECTORY}/05a_final_orthologs_table/final_orthologs_table.csv"
OG_TABLE_ORTHOXML = f"{ALL_OUTPUTS_DIRECTORY}/05a_final_orthologs_table/orthologs.orthoxml"
PHYLETIC_PATTERN = f"{ALL_OUTPUTS_DIRECTORY}/05a_final_orthologs_table/phyletic_pattern.fas"
OG_SIZE_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/05b_groups_sizes_frequency/groups_sizes.png"

GENOME_NUMERIC_REPRESENTATION = f"{ALL_OUTPUTS_DIRECTORY}/06_genome_numeric_representation/core_genome_numeric_representation.txt"

SPECIES_TREE_NEWICK = f"{ALL_OUTPUTS_DIRECTORY}/09_species_phylogeny/final_species_tree.newick"
SPECIES_TREE_PNG = f"{ALL_OUTPUTS_DIRECTORY}/09_species_phylogeny/final_species_tree.png"

OG_TABLE_WITH_CODODN_BIAS = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/final_orthologs_table.csv"
CAI_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/CAI_histogram.png"
GENOMES_CLUSTERS_BY_W_VECTORS = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/Relative_Adaptiveness_scatter_plot.png"
GENOMES_CLUSTERS_BY_W_VECTORS_CSV = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/Relative_Adaptiveness_scatter_plot_clusters.csv"

TITLE_HISTORGRAM_FOR_ORFS = "Open Reading Frames (ORFs)"
TITLE_HISTORGRAM_FOR_GC_CONTENT = "GC Content"
TITLE_HISTOGRAM_GENOME_COMPLETENESS = "Genome Completeness (BUSCO)"
TITLE_HISTOGRAM_ORPHAN_GENES_COUNT = "Orphan Genes Count"

DATA_2_VIEW_IN_HISTOGRAM = {
    # use title as key and path to json as value
    # REMEMBER! json file should be a dict with the genome name as key and a scalar as the value
    TITLE_HISTORGRAM_FOR_ORFS: ORFS_COUNT_PER_GENOME,
    TITLE_HISTORGRAM_FOR_GC_CONTENT: GC_CONTENT_PER_GENOME,
    TITLE_HISTOGRAM_GENOME_COMPLETENESS: GENOME_COMPLETENESS_PER_GENOME,
    TITLE_HISTOGRAM_ORPHAN_GENES_COUNT: ORPHAN_GENES_PER_GENOME
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

source /lsweb/pupko/microbializer/miniconda3/etc/profile.d/conda.sh
conda activate microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

echo "PATH: $PATH"

python "/lsweb/pupko/microbializer/pipeline/main.py" --{args_json_path_key} {args_json_path} --account_name pupkoweb-users --queue_name pupkoweb
echo OKAY > {results_file_path}
'''
