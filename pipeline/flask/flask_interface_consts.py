
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

CORE_PROTEOME = f"{ALL_OUTPUTS_DIRECTORY}/07a_aligned_core_proteome/aligned_core_proteome.fas"
CORE_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/07b_aligned_core_genome/aligned_core_genome.fas"
GENOME_NUMERIC_REPRESENTATION = f"{ALL_OUTPUTS_DIRECTORY}/08_genome_numeric_representation/core_genome_numeric_representation.txt"

SPECIES_TREE_NEWICK = f"{ALL_OUTPUTS_DIRECTORY}/09_species_phylogeny/final_species_tree.newick"
SPECIES_TREE_PNG = f"{ALL_OUTPUTS_DIRECTORY}/09_species_phylogeny/final_species_tree.png"

OG_TABLE_WITH_CODON_BIAS = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/final_orthologs_table.csv"
CAI_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/CAI_histogram.png"
GENOMES_CLUSTERS_BY_W_VECTORS = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/Relative_Adaptiveness_scatter_plot.png"
GENOMES_CLUSTERS_BY_W_VECTORS_CSV = f"{ALL_OUTPUTS_DIRECTORY}/10_codon_bias/Relative_Adaptiveness_scatter_plot_clusters.csv"

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
    "All outputs": {
        "All outputs (zip)": ALL_OUTPUTS_ZIPPED_FORMAT
    },
    "ANI (Average Nucleotide Identity": {
        "ANI table": ANI_CSV,
        "ANI map": ANI_MAP,
    },
    "Genomes statistics": {
        "ORFs count per genome": ORFS_COUNT_PER_GENOME,
        "ORFs count histogram": ORFS_COUNT_HISTOGRAM,
        "GC content per genome": GC_CONTENT_PER_GENOME,
        "GC content histogram": GC_CONTENT_HISTOGRAM,
        "Genome completeness score (BUSCO) per genome": GENOME_COMPLETENESS_PER_GENOME,
        "Genome completeness score (BUSCO) histogram": GENOME_COMPLETENESS_HISTOGRAM,
        "Orphan genes per genome": ORPHAN_GENES_PER_GENOME,
        "Orphan genes histogram": ORPHAN_GENES_HISTOGRAM,
    },
    "Orthologs groups": {
        "Orthologs groups (csv)": OG_TABLE,
        "Orthologs groups (OrthoXML)": OG_TABLE_ORTHOXML,
        "Phyletic pattern": PHYLETIC_PATTERN,
        "Orthologs groups sizes histogram": OG_SIZE_HISTOGRAM,
    },
    "Core genome": {
        "Core genome alignment": CORE_GENOME,
        "Core proteome alignment": CORE_PROTEOME,
        "Genome numeric representation": GENOME_NUMERIC_REPRESENTATION,
    },
    "Species tree": {
        "Species tree (newick)": SPECIES_TREE_PNG,
        "Species tree (png)": SPECIES_TREE_PNG,
    },
    "Codon bias analysis": {
        "Orthologs groups (csv) sorted by codon adaptation index": OG_TABLE_WITH_CODON_BIAS,
        "Codon adaptation index histogram": CAI_HISTOGRAM,
        "Genome clusters by codon usage (png)": GENOMES_CLUSTERS_BY_W_VECTORS,
        "Genome clusters by codon usage (csv)": GENOMES_CLUSTERS_BY_W_VECTORS_CSV
    }
}

# Microbializer processor Job variables
MICROBIALIZER_PROCESSOR_JOB_QUEUE_NAME = 'pupkoweb'
NUBMER_OF_CPUS_MICROBIALIZER_PROCESSOR_JOB = '1'
MICROBIALIZER_PROCESSOR_JOB_PREFIX = 'MC'
MICROBIALIZER_PROCESSOR_RESULTS_FILE_NAME = ALL_OUTPUTS_ZIPPED_FORMAT

MICROBIALIZER_JOB_TEMPLATE = '''#!/bin/bash

sleep {sleep_interval}

echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"

source /lsweb/pupko/microbializer/miniconda3/etc/profile.d/conda.sh
conda activate /lsweb/pupko/microbializer/miniconda3/envs/microbializer
export PATH=$CONDA_PREFIX/bin:$PATH

echo "PATH: $PATH"

python "/lsweb/pupko/microbializer/pipeline/main.py" --{args_json_path_key} {args_json_path} --account_name pupkoweb-users --queue_name pupkoweb
'''

MICROBIALIZER_JOB_HEADER_TEMPLATE = '''
#SBATCH --job-name=microbializer
#SBATCH --account=pupkoweb-users
#SBATCH --partition=pupkoweb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={num_cpus}
#SBATCH --mem=4G
#SBATCH --output={output_file}
#SBATCH --error={error_file}
'''
