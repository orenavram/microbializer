
WEBSERVER_NAME = 'M1CR0B1AL1Z3R'

# Arguments keys to run the pipeline with
CONTIGS_DIR = "contigs_dir"  # can be a path of a directory or a zipped file
IDENTITY_CUTOFF = "identity_cutoff"
E_VALUE_CUTOFF = "e_value_cutoff"
CORE_MINIMAL_PERCENTAGE = "core_minimal_percentage"
BOOTSTRAP = "bootstrap"
OUTGROUP = "outgroup"
FILTER_OUT_PLASMIDS = "filter_out_plasmids"
INPUTS_ARE_ANNOTATED_PROTEOMES = "inputs_are_annotated_proteomes"

# Output files (the paths are relative to the unique folder of the job)
ERROR_FILE_PATH = "error.txt"
ALL_OUTPUTS_ZIPPED_FORMAT = WEBSERVER_NAME + '_{}_outputs'  # need to be formatted with the run_number (which is the job directory name)
