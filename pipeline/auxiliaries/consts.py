#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os.path
from enum import Enum
from pathlib import Path

# ENV = 'wsl'
ENV = 'yair_test'
# ENV = 'yair_prod'
# ENV = 'lsweb'
# ENV = 'c-001'

if ENV == 'yair_test':
    PROJECT_ROOT_DIR = Path('/groups/pupko/yairshimony/microbializer')
elif ENV == 'yair_prod':
    PROJECT_ROOT_DIR = Path('/groups/pupko/yairshimony/microbializer_prod')
elif ENV == 'lsweb':
    PROJECT_ROOT_DIR = Path('/lsweb/pupko/microbializer')
    USER_RESULTS_DIR = PROJECT_ROOT_DIR / 'user_results'
    CLEAN_JOBS_LOGS_DIR = PROJECT_ROOT_DIR / 'clean_jobs_logs'
elif ENV == 'wsl':
    PROJECT_ROOT_DIR = Path('/home/yair/microbializer')
elif ENV == 'c-001':
    PROJECT_ROOT_DIR = Path('/home/ai_center/ai_users/yairshimony/microbializer/')
else:
    raise ValueError(f'Unknown environment: {ENV}')

SRC_DIR = PROJECT_ROOT_DIR / 'pipeline'

if ENV == 'yair_test' or ENV == 'yair_prod':
    CONDA_INSTALLATION_DIR = Path('/groups/pupko/yairshimony/miniconda3')
    CONDA_ENVIRONMENT_DIR = Path('/groups/pupko/yairshimony/miniconda3/envs/microbializer')
    SLURM_ACCOUNT = 'pupko-users'
    SLURM_PARTITION = 'pupko'
elif ENV == 'lsweb':
    CONDA_INSTALLATION_DIR = Path('/lsweb/pupko/microbializer/miniconda3')
    CONDA_ENVIRONMENT_DIR = Path('/lsweb/pupko/microbializer/miniconda3/envs/microbializer')
    SLURM_ACCOUNT = 'pupkoweb-users'
    SLURM_PARTITION = 'pupkoweb'
elif ENV == 'wsl':
    CONDA_INSTALLATION_DIR = Path('/home/yair/miniconda3')
    CONDA_ENVIRONMENT_DIR = Path('/home/yair/miniconda3/envs/microbializer')
    SLURM_ACCOUNT = None
    SLURM_PARTITION = None
elif ENV == 'c-001':
    CONDA_INSTALLATION_DIR = Path('/home/ai_center/ai_users/yairshimony/miniconda')
    CONDA_ENVIRONMENT_DIR = Path('/home/ai_center/ai_users/yairshimony/miniconda/envs/microbializer')
    SLURM_ACCOUNT = 'gpu-research'
    SLURM_PARTITION = 'cpu-killable'
else:
    raise ValueError(f'Unknown environment: {ENV}')

USE_JOB_MANAGER = False if ENV == 'wsl' else True
MAX_PARALLEL_JOBS = 50


# General Job submission consts
Q_SUBMITTER_ADD_SSH_PREFIX = False
if ENV == 'c-001':
    LOGIN_NODE = 'c-001'
else:
    LOGIN_NODE = 'powerslurm-login'

JOB_NAME_ENVIRONMENT_VARIABLE = 'SLURM_JOB_NAME'
JOB_ID_ENVIRONMENT_VARIABLE = 'SLURM_JOB_ID'
JOB_FILES_DEBUG_MODE = False

JOB_WALL_TIME_KEY ='RunTime='
JOB_CPUS_KEY = 'NumCPUs='

DEFAULT_MEMORY_PER_JOB_GB = '4'
MMSEQS_BIG_DATASET_NUM_OF_CORES = 20
MMSEQS_BIG_DATASET_REQUIRED_MEMORY_GB = '64'
PHYLOGENY_NUM_OF_CORES = 20
PHYLOGENY_REQUIRED_MEMORY_GB = '64'
KEGG_NUM_OF_CORES = 20
KEGG_REQUIRED_MEMORY_GB = '64'
CODON_BIAS_NUM_OF_CORES = 20
ANI_NUM_OF_CORES = 20
ANI_REQUIRED_MEMORY_GB = '64'
ORTHOXML_REQUIRED_MEMORY_GB = '16'

INFER_ORTHOGROUPS_JOB_TIME_LIMIT_HOURS = 120
MMSEQS_JOB_TIME_LIMIT_HOURS = 96
PHYLOGENY_JOB_TIME_LIMIT_HOURS = 96

HEGS_ECOLI_FILE_PATH = SRC_DIR / 'data' / 'HEG_ecoli.txt'
KEGG_DATABASE_PATH = SRC_DIR / 'data' / 'kegg' / 'prokaryote_database.hmm'
KEGG_KO_LIST_PATH = SRC_DIR / 'data' / 'kegg' / 'ko_list'

BACTERIA_CORE_GENES_HMM_PROFILES_PATH = SRC_DIR / 'data' / 'busco_hmms'
BACTERIA_CORE_GENES_COUNT = float(len(list(BACTERIA_CORE_GENES_HMM_PROFILES_PATH.iterdir())))
BUSCO_EVAULE_CUTOFF = 10 ** (-2)

MIN_NUMBER_OF_GENOMES_TO_ANALYZE = 2
MAX_NUMBER_OF_GENOMES_TO_ANALYZE = 1060

NUMBER_OF_IQTREE_BOOTSTRAP_ITERATIONS = 1000
NUMBER_OF_RAXML_BOOTSTRAP_ITERATIONS = 100
MAX_NUMBER_OF_CORE_OGS_FOR_PHYLOGENY = 1000
OUTPUT_TSV_OF_ORTHOLOGS_PAIRS = False

DEFAULT_IDENTITY_CUTOFF = 40
DEFAULT_COVERAGE_CUTOFF = 70
DEFAULT_EVALUE_CUTOFF = 0.01
DEFAULT_MMSEQS_SENSITIVITY = 5.7
DEFAULT_CORE_MINIMAL_PERCENTAGE = 100

JOB_DONE_FILE_SUFFIX = '.done'
CHECK_JOB_DONE_INTERVAL_SECONDS = 10

class SimilarityScore(Enum):
    BITS = 1
    EVALUE = 2  # Currently not supported


SIMILARITY_SCORE_CRITERION = SimilarityScore.BITS

BLAST_OUTPUT_HEADER = ['query', 'subject', 'identity_percent', 'alignment_length', 'mismatches', 'gap_openings',
                        'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bit_score']
MMSEQS_OUTPUT_FORMAT = 'query,target,bits'
MMSEQS_OUTPUT_HEADER = MMSEQS_OUTPUT_FORMAT.split(',')
MMSEQS_OUTPUT_COLUMNS_TYPES = {'query': str, 'target': str, 'bits': int}
HMMSEARCH_OUTPUT_HEADER = ['target_name', 'target_accession', 'query_name', 'query_accession', 'full_e_value',
                           'full_score', 'full_bias', 'domain_e_Value', 'domain_score', 'domain_bias', 'exp', 'reg',
                           'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']

# logging consts
LOG_MESSAGE_FORMAT = '%(asctime)s %(levelname)s %(filename)s:%(lineno)d %(message)s'

# outputs directories map
OUTPUTS_DIRECTORIES_MAP = {
    '11_1_ani': '01_ani',
    'orfs_sequences': '02a_orfs',
    '02_2_orfs_plots': '02b_orfs_plots',
    'orfs_translated': '02c_translated_orfs',
    '03_genomes_completeness': '03_genomes_completeness',
    '05_9_orphan_genes': '04_orphan_genes',    # full_orthogroups_inference
    '05_10_orthogroups_final': '05a_orthogroups',   # full_orthogroups_inference
    '06_12_orphan_genes_from_orthogroups': '04_orphan_genes',    # approximate_orthogroups_inference
    '06_13_orthogroups_final': '05a_orthogroups',    # approximate_orthogroups_inference
    '07_1_orthogroups_variations': '05a_orthogroups',
    '07_2_orthogroups_sizes': '05b_orthogroups_sizes',
    'orthogroups_dna': '06a_orthogroups_dna',
    'orthogroups_aa': '06b_orthogroups_aa',
    'orthogroups_aa_msa': '06c_orthogroups_aa_msa',
    'orthogroups_induced_dna_msa': '06d_orthogroups_induced_dna_msa_by_aa_msa',
    '09_1_aligned_core_proteome': '07a_aligned_core_proteome',
    '09_2_aligned_core_genome': '07b_aligned_core_genome',
    '10_genome_numeric_representation': '08_genome_numeric_representation',
    '11_2_species_phylogeny': '09_species_phylogeny',
    '12_1_codon_bias': '10_codon_bias',
    '12_3_orthogroups_annotations': '05a_orthogroups',
}

# steps for progress bar
FULL_STEPS_NAMES_FOR_PROGRESS_BAR = [
    'Validate input files',
    'Filter out plasmids',
    'Predict and translate ORFs',
    'Calculate genomes completeness',
    'Infer orthogroups',
    'Find orphan genes',
    'Prepare orthogroups fasta files',
    'Infer core genome and proteome',
    'Calculate genomes numeric representation',
    'Reconstruct species phylogeny',
    'Calculate ANI (Average Nucleotide Identity)',
    'Analyze orthogroups codon bias',
    'Annotate orthogroups with KO terms',
    'Finalize results',
]

ONLY_CALC_OGS_TABLE_STEPS_NAMES_FOR_PROGRESS_BAR = [
    'Validate input files',
    'Filter out plasmids',
    'Predict and translate ORFs',
    'Infer orthogroups',
    'Find orphan genes',
    'Finalize results',
]
