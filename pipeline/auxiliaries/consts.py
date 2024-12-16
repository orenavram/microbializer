#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os.path
from enum import Enum
import numpy as np

KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR = True
LOG_IN_SEPARATE_FILES = True

# ENV = 'wsl'
ENV = 'yair_test'
# ENV = 'yair_prod'
# ENV = 'lsweb'
# ENV = 'c-001'

if ENV == 'yair_test':
    PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer'
elif ENV == 'yair_prod':
    PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer_prod'
elif ENV == 'lsweb':
    PROJECT_ROOT_DIR = '/lsweb/pupko/microbializer'
elif ENV == 'wsl':
    PROJECT_ROOT_DIR = '/home/yair/microbializer'
elif ENV == 'c-001':
    PROJECT_ROOT_DIR = '/home/ai_center/ai_users/yairshimony/microbializer/'
else:
    raise ValueError(f'Unknown environment: {ENV}')

SRC_DIR = os.path.join(PROJECT_ROOT_DIR, 'pipeline')

if ENV == 'yair_test' or ENV == 'yair_prod':
    CONDA_INSTALLATION_DIR = r'/groups/pupko/yairshimony/miniconda3'
    CONDA_ENVIRONMENT_DIR = r'/groups/pupko/yairshimony/miniconda3/envs/microbializer'
elif ENV == 'lsweb':
    CONDA_INSTALLATION_DIR = r'/lsweb/pupko/microbializer/miniconda3'
    CONDA_ENVIRONMENT_DIR = r'/lsweb/pupko/microbializer/miniconda3/envs/microbializer'
elif ENV == 'wsl':
    CONDA_INSTALLATION_DIR = r'/home/yair/miniconda3'
    CONDA_ENVIRONMENT_DIR = r'/home/yair/miniconda3/envs/microbializer'
elif ENV == 'c-001':
    CONDA_INSTALLATION_DIR = r'/home/ai_center/ai_users/yairshimony/miniconda'
    CONDA_ENVIRONMENT_DIR = r'/home/ai_center/ai_users/yairshimony/miniconda/envs/microbializer'
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

MMSEQS_CLUSTER_MIN_SEQ_ID = 5
MMSEQS_CLUSTER_MIN_COVERAGE = 10
JOB_WALL_TIME_KEY ='RunTime='
JOB_CPUS_KEY = 'NumCPUs='

DEFAULT_MEMORY_PER_JOB_GB = '4'
MMSEQS_BIG_DATASET_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1
MMSEQS_BIG_DATASET_REQUIRED_MEMORY_GB = '64'
MMSEQS_CLUSTER_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1
MMSEQS_CLUSTER_REQUIRED_MEMORY_GB = '64'
PHYLOGENY_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1
PHYLOGENY_REQUIRED_MEMORY_GB = '64'
KEGG_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1
KEGG_REQUIRED_MEMORY_GB = '64'
CODON_BIAS_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1
ANI_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1
ANI_REQUIRED_MEMORY_GB = '64'
MCL_NUM_OF_CORES = 20 if USE_JOB_MANAGER else 1

INFER_ORTHOGROUPS_JOB_TIME_LIMIT_HOURS = 144
MMSEQS_JOB_TIME_LIMIT_HOURS = 72


# Slurm consts
DEFAULT_SLURM_ACCOUNT = 'pupko-users'
DEFAULT_SLURM_PARTITION = 'pupko'

HEGS_ECOLI_FILE_PATH = os.path.join(SRC_DIR, 'data', 'HEG_ecoli.txt')
BACTERIA_CORE_GENES_HMM_PROFILES_PATH = os.path.join(SRC_DIR, 'data', 'busco_hmms')
KEGG_DATABASE_PATH = os.path.join(SRC_DIR, 'data', 'kegg', 'prokaryote_database.hmm')
KEGG_KO_LIST_PATH = os.path.join(SRC_DIR, 'data', 'kegg', 'ko_list')
MAX_NUMBER_OF_GENOMES_TO_ANALYZE = 350
NUMBER_OF_IQTREE_BOOTSTRAP_ITERATIONS = 1000
NUMBER_OF_RAXML_BOOTSTRAP_ITERATIONS = 100
NAX_NUMBER_OF_CORE_OGS_FOR_PHYLOGENY = 1000
OUTPUT_TSV_OF_ORTHOLOGS_PAIRS = False


class SimilarityScore(Enum):
    BITS = 1
    EVALUE = 2


SIMILARITY_SCORE_CRITERION = SimilarityScore.BITS

BLAST_OUTPUT_HEADER = ['query', 'subject', 'identity_percent', 'alignment_length', 'mismatches', 'gap_openings',
                        'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bit_score']
MMSEQS_OUTPUT_FORMAT = 'query,target,bits'
MMSEQS_OUTPUT_HEADER = MMSEQS_OUTPUT_FORMAT.split(',')
MMSEQS_OUTPUT_COLUMNS_TYPES = {'query': str, 'target': str, 'bits': np.float16}
HMMSEARCH_OUTPUT_HEADER = ['target_name', 'target_accession', 'query_name', 'query_accession', 'full_e_value',
                           'full_score', 'full_bias', 'domain_e_Value', 'domain_score', 'domain_bias', 'exp', 'reg',
                           'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']

# logging consts
LOG_MESSAGE_FORMAT = '%(asctime)s %(levelname)s %(filename)s:%(lineno)d %(message)s'

# outputs directories map
OUTPUTS_DIRECTORIES_MAP = {
    '01_ani': '01_ani',
    '02_1_orfs': '02a_orfs',
    '02_3_orfs_plots': '02b_orfs_plots',
    '02_4_translated_orfs': '02c_translated_orfs',
    '03_genomes_completeness': '03_genomes_completeness',
    '06_orphan_genes': '04_orphan_genes',
    '07_1_orthogroups': '05a_orthogroups',
    '07_2_orthogroups_sizes': '05b_orthogroups_sizes',
    '08_1_orthogroups_dna': '06a_orthogroups_dna',
    '08_2_orthogroups_aa': '06b_orthogroups_aa',
    '08_3_orthogroups_aa_msa': '06c_orthogroups_aa_msa',
    '08_4_orthogroups_induced_dna_msa_by_aa_msa': '06d_orthogroups_induced_dna_msa_by_aa_msa',
    '09_1_aligned_core_proteome': '07a_aligned_core_proteome',
    '09_2_aligned_core_genome': '07b_aligned_core_genome',
    '10_genome_numeric_representation': '08_genome_numeric_representation',
    '11_species_phylogeny': '09_species_phylogeny',
    '12_codon_bias': '10_codon_bias',
    'orthogroups_annotated.csv': '05a_orthogroups',  # copy the annotated OG table (with KEGG+codon bias) to the same directory as the simple OG table
}

# steps for progress bar
FULL_STEPS_NAMES_FOR_PROGRESS_BAR = [
    'Validate input files',
    'Filter out plasmids',
    'Calculate ANI (Average Nucleotide Identity)',
    'Predict and translate ORFs',
    'Calculate genomes completeness',
    'Infer orthogroups',
    'Find orphan genes',
    'Prepare orthogroups fasta files',
    'Infer core genome and proteome',
    'Calculate genomes numeric representation',
    'Reconstruct species phylogeny',
    'Analyze codon bias',
    'Annotate orthogroups with KEGG Orthology (KO) terms',
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


# general vars
CSV_DELIMITER = ','
