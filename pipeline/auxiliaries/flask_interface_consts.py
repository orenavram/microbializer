
TEST = True
WEBSERVER_NAME = 'M1CR0B1AL1Z3R'
QUEUE_FOR_JOBS = 'power-pupko'
PROJECT_ROOT_DIR = '/groups/pupko/yairshimony/microbializer' if TEST else '/bioseq/microbializer'
OWNER_EMAIL = 'yairshimony@mail.tau.ac.il' if TEST else 'orenavram@gmail.com'

# Arguments keys to run the pipeline with
CONTIGS_DIR = "contigs_dir"  # can be a path of a directory or a zipped file
IDENTITY_CUTOFF = "identity_cutoff"
E_VALUE_CUTOFF = "e_value_cutoff"
CORE_MINIMAL_PERCENTAGE = "core_minimal_percentage"
BOOTSTRAP = "bootstrap"
OUTGROUP = "outgroup"
FILTER_OUT_PLASMIDS = "filter_out_plasmids"
INPUTS_ARE_ANNOTATED_PROTEOMES = "inputs_are_annotated_proteomes"

# Error description file (the path is relative to the unique folder of the job)
ERROR_FILE_PATH = "error.txt"

# Output files (the paths are relative to the unique folder of the job).
# Need to be formatted with the run_number (which is the job directory name)
ALL_OUTPUTS_DIRECTORY = WEBSERVER_NAME + "_{}_outputs"
ALL_OUTPUTS_ZIPPED_FORMAT = WEBSERVER_NAME + "_{}_outputs.zip"
ORFS_COUNT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_counts.json"
ORFS_COUNT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_counts.png"
GC_CONTENT_PER_GENOME = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_gc_contents.json"
GC_CONTENT_HISTOGRAM = f"{ALL_OUTPUTS_DIRECTORY}/20_orfs_plot/orfs_gc_contents.png"
