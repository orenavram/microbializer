#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os.path
from enum import Enum

KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR = True
USE_CONDA = True
IGNORE_HTML = True
SEND_MAILS = False
LOG_IN_SEPARATE_FILES = True

# ENV = 'yair_test'
ENV = 'yair_prod'
# ENV = 'lsweb'

if ENV == 'yair_test':
    PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer'
elif ENV == 'yair_prod':
    PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer_prod'
elif ENV == 'lsweb':
    PROJECT_ROOT_DIR = '/lsweb/pupko/microbializer'
else:
    raise ValueError(f'Unknown environment: {ENV}')

SRC_DIR = os.path.join(PROJECT_ROOT_DIR, 'pipeline')

if ENV == 'yair_test' or ENV == 'yair_prod':
    CONDA_INSTALLATION_DIR = r'/groups/pupko/yairshimony/miniconda3'
    CONDA_ENVIRONMENT_DIR = r'/groups/pupko/yairshimony/miniconda3/envs/microbializer'
elif ENV == 'lsweb':
    CONDA_INSTALLATION_DIR = r'/lsweb/pupko/microbializer/miniconda3'
    CONDA_ENVIRONMENT_DIR = r'/lsweb/pupko/microbializer/miniconda3/envs/microbializer'
else:
    raise ValueError(f'Unknown environment: {ENV}')

OWNER_EMAIL = 'yairshimony@mail.tau.ac.il'


# General Job submission consts
Q_SUBMITTER_ADD_SSH_PREFIX = False

JOB_NAME_ENVIRONMENT_VARIABLE = 'SLURM_JOB_NAME'
JOB_ID_ENVIRONMENT_VARIABLE = 'SLURM_JOB_ID'
JOB_FILES_DEBUG_MODE = False
PHYLOGENY_NUM_OF_CORES = 20
CODON_BIAS_NUM_OF_CORES = 20
CLUSTER_PROTEOMES_NUM_OF_CORES = 20
MMSEQS_CLUSTER_MIN_SEQ_ID = 5
MMSEQS_CLUSTER_MIN_COVERAGE = 10
JOB_WALL_TIME_KEY ='RunTime='
# mmseqs and fastANI commands work only on machines with enough memory. we solve this either by navigating to a
# specific queue or by restrict the compute-nodes with memory threshold.
MMSEQS_REQUIRED_MEMORY_GB = '120'
ANI_REQUIRED_MEMORY_GB = '120'
PHYLOGENY_REQUIRED_MEMORY_GB = '120'

# Slurm consts
DEFAULT_SLURM_ACCOUNT = 'pupko-users'
DEFAULT_SLURM_PARTITION = 'pupko'

HEGS_ECOLI_FILE_PATH = os.path.join(SRC_DIR, 'data', 'HEG_ecoli.txt')
BACTERIA_CORE_GENES_HMM_PROFILES_PATH = os.path.join(PROJECT_ROOT_DIR, 'pipeline', 'data', 'busco_hmms')
MAX_NUMBER_OF_GENOMES_TO_ANALYZE = 350
NUMBER_OF_IQTREE_BOOTSTRAP_ITERATIONS = 1000
NUMBER_OF_RAXML_BOOTSTRAP_ITERATIONS = 100
OUTPUT_TSV_OF_ORTHOLOGS_PAIRS = False


class SimilarityScore(Enum):
    BITS = 1
    EVALUE = 2


SIMILARITY_SCORE_CRITERION = SimilarityScore.BITS

BLAST_OUTPUT_HEADER = ['query', 'subject', 'identity_percent', 'alignment_length', 'mismatches', 'gap_openings',
                        'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bit_score']
MMSEQS_OUTPUT_FORMAT = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,qcov,tstart,tend,tlen,tcov,evalue,bits'
MMSEQS_OUTPUT_HEADER = MMSEQS_OUTPUT_FORMAT.split(',')

# logging consts
LOG_MESSAGE_FORMAT = '%(asctime)s %(levelname)s %(filename)s:%(lineno)d %(message)s'

# outputs directories map
OUTPUTS_DIRECTORIES_MAP = {
    '01_ANI': '01_ANI',
    '02a_orfs': '02a_orfs',
    '02c_orfs_plots': '02b_orfs_plots',
    '02d_translated_orfs': '02c_translated_orfs',
    '03_genomes_completeness': '03_genomes_completeness',
    '06_orphan_genes': '04_orphan_genes',
    '07a_final_orthologs_table': '05a_final_orthologs_table',
    '07b_groups_sizes_frequency': '05b_groups_sizes_frequency',
    '08a_orthologs_groups_dna': '06a_orthologs_groups_dna',
    '08b_orthologs_groups_aa': '06b_orthologs_groups_aa',
    '08c_orthologs_groups_aa_msa': '06c_orthologs_groups_aa_msa',
    '08d_orthologs_groups_induced_dna_msa_by_aa_msa': '06d_orthologs_groups_induced_dna_msa_by_aa_msa',
    '09a_aligned_core_proteome': '07a_aligned_core_proteome',
    '09b_aligned_core_genome': '07b_aligned_core_genome',
    '10_genome_numeric_representation': '08_genome_numeric_representation',
    '11_species_phylogeny': '09_species_phylogeny',
    '12_codon_bias': '10_codon_bias',
}

# steps for progress bar
FULL_STEPS_NAMES_FOR_PROGRESS_BAR = [
    'Validate input files',
    'Filter out plasmids',
    'Calculate ANI (Average Nucleotide Identity)',
    'Predict and translate ORFs',
    'Calculate genomes completeness',
    'Calculate blast scores between all protein sequences',
    'Cluster protein sequences',
    'Find orphan genes',
    'Construct orthogroups',
    'Prepare orthogroups fasta files',
    'Infer core genome',
    'Calculate genomes numeric representation',
    'Reconstruct species phylogeny',
    'Analyze codon bias',
    'Finalize results',
]

ONLY_CALC_OGS_TABLE_STEPS_NAMES_FOR_PROGRESS_BAR = [
    'Validate input files',
    'Filter out plasmids',
    'Predict and translate ORFs',
    'Calculate blast scores between all protein sequences',
    'Cluster protein sequences',
    'Find orphan genes',
    'Construct orthogroups',
    'Finalize results',
]


# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

# general vars
CSV_DELIMITER = ','

# relevant modules
GCC = 'gcc/gcc-6.2.0'
MCL = 'MCL-edge/mcl-14-137'
MAFFT = 'mafft/mafft-7.407'
RAXML = 'raXML'
PRODIGAL = 'prodigal/prodigal-2.6.3'
MMSEQS = 'MMseqs2/June2020'

