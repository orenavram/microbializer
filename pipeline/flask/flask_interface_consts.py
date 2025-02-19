from pathlib import Path

WEBSERVER_NAME = 'M1CR0B1AL1Z3R'
WEBSERVER_PROJECT_ROOT_DIR = '/lsweb/pupko/microbializer'
ADDITIONAL_OWNER_EMAILS = ['edodotan@mail.tau.ac.il']

# Arguments keys to run the pipeline with
ARGS_JSON_PATH_KEY = "args_json_path"
JOB_PARAMETERS_FILE_NAME = "input_parameters.json"
RUN_DIR = "run_dir"

# Input parameters from users
JOB_NAME = "job_name"
EMAIL = "email"
IDENTITY_CUTOFF = "identity_cutoff"
E_VALUE_CUTOFF = "e_value_cutoff"
CORE_MINIMAL_PERCENTAGE = "core_minimal_percentage"
COVERAGE_CUTOFF = "coverage_cutoff"
BOOTSTRAP = "bootstrap"
OUTGROUP = "outgroup"
FILTER_OUT_PLASMIDS = "filter_out_plasmids"
ADD_ORPHAN_GENES_TO_OGS = "add_orphan_genes_to_ogs"
INPUT_FASTA_TYPE = "inputs_fasta_type"

# Input file
INPUTS_GENOMES_ZIPPED = "genomes.zip"

# Output files (the paths are relative to the unique folder of the job)
ALL_OUTPUTS_DIRECTORY = Path(WEBSERVER_NAME + "_outputs")
ALL_OUTPUTS_ZIPPED = WEBSERVER_NAME + "_outputs.zip"

# Error description file (the path is relative to the unique folder of the job)
ERROR_FILE_NAME = "error.txt"
ERROR_FILE_PATH = ALL_OUTPUTS_DIRECTORY / ERROR_FILE_NAME

# Whether to send email when job finished from pipeline or flask
SEND_EMAIL_WHEN_JOB_FINISHED_FROM_PIPELINE = True  # If True, flask won't send emails. If False, flask will send emails.
# Whether to clean old jobs directories from pipeline or flask
CLEAN_OLD_JOBS_DIRECTORIES_FROM_PIPELINE = True  # If True, flask won't clean old jobs directories. If False, flask will clean old jobs directories.

# Progress bar file (the path is relative to the unique folder of the job)
PROGRESSBAR_FILE_NAME = "progressbar.csv"

ANI_CSV = ALL_OUTPUTS_DIRECTORY / '01_ani' / 'ani_pairwise_values.csv'
ANI_MAP = ALL_OUTPUTS_DIRECTORY / '01_ani' / 'ani_map.png'

ORFS_COUNT_PER_GENOME = ALL_OUTPUTS_DIRECTORY / '02b_orfs_plots' / 'orfs_counts.json'
ORFS_COUNT_HISTOGRAM = ALL_OUTPUTS_DIRECTORY / '02b_orfs_plots' / 'orfs_counts.png'
GC_CONTENT_PER_GENOME = ALL_OUTPUTS_DIRECTORY / '02b_orfs_plots' / 'orfs_gc_content.json'
GC_CONTENT_HISTOGRAM = ALL_OUTPUTS_DIRECTORY / '02b_orfs_plots' / 'orfs_gc_content.png'
GENOME_COMPLETENESS_PER_GENOME = ALL_OUTPUTS_DIRECTORY / '03_genomes_completeness' / 'genomes_completeness.json'
GENOME_COMPLETENESS_HISTOGRAM = ALL_OUTPUTS_DIRECTORY / '03_genomes_completeness' / 'genomes_completeness.png'
ORPHAN_GENES_PER_GENOME = ALL_OUTPUTS_DIRECTORY / '04_orphan_genes' / 'orphan_genes_count.json'
ORPHAN_GENES_HISTOGRAM = ALL_OUTPUTS_DIRECTORY / '04_orphan_genes' / 'orphan_genes_count.png'

OG_TABLE = ALL_OUTPUTS_DIRECTORY / '05a_orthogroups' / 'orthogroups.csv'
OG_TABLE_ANNOTATED = ALL_OUTPUTS_DIRECTORY / '05a_orthogroups' / 'orthogroups_annotated.csv'
OG_TABLE_ORTHOXML = ALL_OUTPUTS_DIRECTORY / '05a_orthogroups' / 'orthogroups.orthoxml'
PHYLETIC_PATTERN = ALL_OUTPUTS_DIRECTORY / '05a_orthogroups' / 'phyletic_pattern.fas'
OG_SIZE_HISTOGRAM = ALL_OUTPUTS_DIRECTORY / '05b_orthogroups_sizes' / 'groups_sizes.png'

CORE_PROTEOME = ALL_OUTPUTS_DIRECTORY / '07a_aligned_core_proteome' / 'aligned_core_proteome.fas'
CORE_GENOME = ALL_OUTPUTS_DIRECTORY / '07b_aligned_core_genome' / 'aligned_core_genome.fas'
GENOME_NUMERIC_REPRESENTATION = ALL_OUTPUTS_DIRECTORY / '08_genome_numeric_representation' / 'core_genome_numeric_representation.txt'

SPECIES_TREE_NEWICK = ALL_OUTPUTS_DIRECTORY / '09_species_phylogeny' / 'final_species_tree.newick'
SPECIES_TREE_PNG = ALL_OUTPUTS_DIRECTORY / '09_species_phylogeny' / 'final_species_tree.png'

CAI_HISTOGRAM = ALL_OUTPUTS_DIRECTORY / '10_codon_bias' / 'CAI_histogram.png'
W_VECTORS = ALL_OUTPUTS_DIRECTORY / '10_codon_bias' / 'W_vectors.csv'
GENOMES_CLUSTERS_BY_W_VECTORS = ALL_OUTPUTS_DIRECTORY / '10_codon_bias' / 'Relative_Adaptiveness_scatter_plot.png'
GENOMES_CLUSTERS_BY_W_VECTORS_CSV = ALL_OUTPUTS_DIRECTORY / '10_codon_bias' / 'Relative_Adaptiveness_scatter_plot_clusters.csv'

TITLE_HISTOGRAM_FOR_ORFS = "Open Reading Frames (ORFs)"
TITLE_HISTOGRAM_FOR_GC_CONTENT = "GC Content %"
TITLE_HISTOGRAM_GENOME_COMPLETENESS = "Genome Completeness (BUSCO)"
TITLE_HISTOGRAM_ORPHAN_GENES_COUNT = "Orphan Genes Count"

DATA_2_VIEW_IN_HISTOGRAM = {
    # use title as key and path to json as value
    # REMEMBER! json file should be a dict with the genome name as key and a scalar as the value
    TITLE_HISTOGRAM_FOR_ORFS: ORFS_COUNT_PER_GENOME,
    TITLE_HISTOGRAM_FOR_GC_CONTENT: GC_CONTENT_PER_GENOME,
    TITLE_HISTOGRAM_GENOME_COMPLETENESS: GENOME_COMPLETENESS_PER_GENOME,
    TITLE_HISTOGRAM_ORPHAN_GENES_COUNT: ORPHAN_GENES_PER_GENOME
}

PATHS_TO_DOWNLOAD = {
    "General": {
        "Inputs_(zip)": (INPUTS_GENOMES_ZIPPED, "The input genomes in a zip file"),
        "All_outputs_(zip)": (ALL_OUTPUTS_ZIPPED, "All the outputs in a zip file"),
    },
    "ANI (Average Nucleotide Identity)": {
        "ANI_table": (ANI_CSV, "Can also be found in 01_ani directory in the 'All outputs (zip)' file"),
        "ANI_map": (ANI_MAP, "Can also be found in 01_ani directory in the 'All outputs (zip)' file"),
    },
    "Genomes statistics": {
        "ORFs_count_per_genome": (ORFS_COUNT_PER_GENOME, "Can also be found in 02b_orfs_plots directory in the 'All outputs (zip)' file"),
        "ORFs_count_histogram": (ORFS_COUNT_HISTOGRAM, "Can also be found in 02b_orfs_plots directory in the 'All outputs (zip)' file"),
        "GC_content_per_genome": (GC_CONTENT_PER_GENOME, "Can also be found in 02b_orfs_plots directory in the 'All outputs (zip)' file"),
        "GC_content_histogram": (GC_CONTENT_HISTOGRAM, "Can also be found in 02b_orfs_plots directory in the 'All outputs (zip)' file"),
        "Genome_completeness_score_(BUSCO)_per_genome": (GENOME_COMPLETENESS_PER_GENOME, "Can also be found in 03_genomes_completeness directory in the 'All outputs (zip)' file"),
        "Genome_completeness_score_(BUSCO)_histogram": (GENOME_COMPLETENESS_HISTOGRAM, "Can also be found in 03_genomes_completeness directory in the 'All outputs (zip)' file"),
        "Orphan_genes_per_genome": (ORPHAN_GENES_PER_GENOME, "Can also be found in 04_orphan_genes directory in the 'All outputs (zip)' file"),
        "Orphan_genes_histogram": (ORPHAN_GENES_HISTOGRAM, "Can also be found in 04_orphan_genes directory in the 'All outputs (zip)' file"),
    },
    "Orthogroups": {
        "Orthogroups_(csv)": (OG_TABLE, "Can also be found in 05a_orthogroups directory in the 'All outputs (zip)' file"),
        "Orthogroups_annotated_(csv)": (OG_TABLE_ANNOTATED, "Can also be found in 05a_orthogroups directory in the 'All outputs (zip)' file"),
        "Orthogroups_(OrthoXML)": (OG_TABLE_ORTHOXML, "Can also be found in 05a_orthogroups directory in the 'All outputs (zip)' file"),
        "Phyletic_pattern": (PHYLETIC_PATTERN, "Can also be found in 05a_orthogroups directory in the 'All outputs (zip)' file"),
        "Orthogroups_sizes_histogram": (OG_SIZE_HISTOGRAM, "Can also be found in 05b_orthogroups_sizes directory in the 'All outputs (zip)' file"),
    },
    "Core genome": {
        "Core_proteome_alignment": (CORE_GENOME, "Can also be found in 07a_aligned_core_proteome directory in the 'All outputs (zip)' file"),
        "Core_genome_alignment": (CORE_PROTEOME, "Can also be found in 07b_aligned_core_genome directory in the 'All outputs (zip)' file"),
        "Genome_numeric_representation": (GENOME_NUMERIC_REPRESENTATION, "Can also be found in 08_genome_numeric_representation directory in the 'All outputs (zip)' file"),
    },
    "Species tree": {
        "Species_tree_(newick)": (SPECIES_TREE_NEWICK, "Can also be found in 09_species_phylogeny directory in the 'All outputs (zip)' file"),
        "Species_tree_(png)": (SPECIES_TREE_PNG, "Can also be found in 09_species_phylogeny directory in the 'All outputs (zip)' file"),
    },
    "Codon bias analysis": {
        "Codon_adaptation_index_histogram": (CAI_HISTOGRAM, "Can also be found in 10_codon_bias directory in the 'All outputs (zip)' file"),
        "W_vectors_(relative_adaptiveness)": (W_VECTORS, "Can also be found in 10_codon_bias directory in the 'All outputs (zip)' file"),
        "Genome_clusters_by_relative_adaptiveness_(png)": (GENOMES_CLUSTERS_BY_W_VECTORS, "Can also be found in 10_codon_bias directory in the 'All outputs (zip)' file"),
        "Genome_clusters_by_relative_adaptiveness_(csv)": (GENOMES_CLUSTERS_BY_W_VECTORS_CSV, "Can also be found in 10_codon_bias directory in the 'All outputs (zip)' file"),
    }
}

# Microbializer processor Job variables
MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME = 'pupkoweb'
MICROBIALIZER_PROCESSOR_JOB_ACCOUNT_NAME = 'pupkoweb-users'
NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB = '1'
MICROBIALIZER_MAIN_JOB_MEMORY = '8192'  # in MB
MICROBIALIZER_MAIN_JOB_TIME_LIMIT_IN_HOURS = 168  # 7 days (the max time of pupkoweb)
MICROBIALIZER_PROCESSOR_JOB_PREFIX = 'MC'
MICROBIALIZER_PROCESSOR_RESULTS_FILE_NAME = ALL_OUTPUTS_ZIPPED
INTERVAL_BETWEEN_LISTENER_SAMPLES = 120  # in seconds

MICROBIALIZER_JOB_TEMPLATE = f'''#!/bin/bash

sleep 30

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Memory per node: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "Job name: $SLURM_JOB_NAME"

export HOME={WEBSERVER_PROJECT_ROOT_DIR}
source {WEBSERVER_PROJECT_ROOT_DIR}/miniconda3/etc/profile.d/conda.sh
conda activate {WEBSERVER_PROJECT_ROOT_DIR}/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

echo "PATH: $PATH"

python "{WEBSERVER_PROJECT_ROOT_DIR}/pipeline/main.py" --{ARGS_JSON_PATH_KEY} {{args_json_path}} --account_name {MICROBIALIZER_PROCESSOR_JOB_ACCOUNT_NAME} --queue_name {MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME}
'''

MICROBIALIZER_JOB_HEADER_TEMPLATE = f'''
#SBATCH --job-name=microbializer
#SBATCH --account={MICROBIALIZER_PROCESSOR_JOB_ACCOUNT_NAME}
#SBATCH --partition={MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB}
#SBATCH --mem={MICROBIALIZER_MAIN_JOB_MEMORY}
#SBATCH --time={MICROBIALIZER_MAIN_JOB_TIME_LIMIT_IN_HOURS}:00:00
#SBATCH --output={{output_file}}
#SBATCH --error={{error_file}}
'''
