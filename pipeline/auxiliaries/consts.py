#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os.path
from enum import Enum

KEEP_OUTPUTS_IN_INTERMEDIATE_RESULTS_DIR = True
USE_TEST_ROOT_DIR = True
USE_CONDA = True
IGNORE_HTML = True
LOG_IN_SEPARATE_FILES = True
PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer' if USE_TEST_ROOT_DIR else \
    '/groups/pupko/yairshimony/microbializer_prod'
CONDA_INSTALLATION_DIR = r'/groups/pupko/yairshimony/miniconda3'
CONDA_ENVIRONMENT_DIR = r'/groups/pupko/yairshimony/miniconda3/envs/microbializer'
OWNER_EMAIL = 'yairshimony@mail.tau.ac.il'


# General Job submission consts
Q_SUBMITTER_ADD_SSH_PREFIX = False
PBS = False  # if False, assume slurm
JOB_NAME_ENVIRONMENT_VARIABLE = 'PBS_JOBNAME' if PBS else 'SLURM_JOB_NAME'
JOB_ID_ENVIRONMENT_VARIABLE = 'PBS_JOBID' if PBS else 'SLURM_JOB_ID'
JOB_FILES_DEBUG_MODE = False
PHYLOGENY_NUM_OF_CORES = 20
CODON_BIAS_NUM_OF_CORES = 20
# mmseqs and fastANI commands work only on machines with enough memory. we solve this either by navigating to a
# specific queue or by restrict the compute-nodes with memory threshold.
MMSEQS_REQUIRED_MEMORY_GB = '120'
ANI_REQUIRED_MEMORY_GB = '120'
USE_DIFFERENT_QUEUE_FOR_MMSEQS = False
QUEUE_FOR_MMSEQS_COMMANDS = 'pupkolab'

# PBS consts
PBS_QUEUE = 'power-pupko'

# Slurm consts
SLURM_ACCOUNT = 'power-general-users'
SLURM_PARTITION = 'power-general'


HEGS_ECOLI_FILE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'HEG_ecoli.txt')
BACTERIA_CORE_GENES_HMM_PROFILES_PATH = '/groups/pupko/naamawagner/Microbializer/Busco/hmms'
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

