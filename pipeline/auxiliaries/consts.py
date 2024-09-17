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
# ENV = 'yair_prod'
ENV = 'lsweb'

if ENV == 'yair_test':
    PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer'
elif ENV == 'yair_prod':
    PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer_prod'
elif ENV == 'lsweb':
    PROJECT_ROOT_DIR = '/lsweb/pupko/microbializer'
else:
    raise ValueError(f'Unknown environment: {ENV}')

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
PBS = False  # if False, assume slurm
JOB_NAME_ENVIRONMENT_VARIABLE = 'PBS_JOBNAME' if PBS else 'SLURM_JOB_NAME'
JOB_ID_ENVIRONMENT_VARIABLE = 'PBS_JOBID' if PBS else 'SLURM_JOB_ID'
JOB_FILES_DEBUG_MODE = False
PHYLOGENY_NUM_OF_CORES = 20
CODON_BIAS_NUM_OF_CORES = 20
JOB_CPU_TIME_KEY = 'resources_used.cput = ' if PBS else '' # Only in PBS I found a way to the get the job's cpu runtime from within the job (in the compute node)
JOB_WALL_TIME_KEY = 'resources_used.walltime = ' if PBS else 'RunTime='
# mmseqs and fastANI commands work only on machines with enough memory. we solve this either by navigating to a
# specific queue or by restrict the compute-nodes with memory threshold.
MMSEQS_REQUIRED_MEMORY_GB = '120'
ANI_REQUIRED_MEMORY_GB = '120'
USE_DIFFERENT_QUEUE_FOR_MMSEQS = False
QUEUE_FOR_MMSEQS_COMMANDS = 'pupkolab'

# PBS consts
DEFAULT_PBS_QUEUE = 'power-pupko'

# Slurm consts
DEFAULT_SLURM_ACCOUNT = 'power-general-users'
DEFAULT_SLURM_PARTITION = 'power-general'

HEGS_ECOLI_FILE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'HEG_ecoli.txt')
BACTERIA_CORE_GENES_HMM_PROFILES_PATH = os.path.join(PROJECT_ROOT_DIR, 'pipeline', 'data', 'busco_hmms')
MAX_NUMBER_OF_GENOMES_TO_ANALYZE = 350
NUMBER_OF_IQTREE_BOOTSTRAP_ITERATIONS = 1000
NUMBER_OF_RAXML_BOOTSTRAP_ITERATIONS = 100


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

