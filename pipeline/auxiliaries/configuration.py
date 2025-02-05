from dataclasses import dataclass, asdict
import argparse
import json
import os
import logging
import sys
from pathlib import Path
import pandas as pd

from . import consts
from .pipeline_auxiliaries import str_to_bool

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parent))

from flask.flask_interface_consts import (WEBSERVER_NAME, ERROR_FILE_NAME, PROGRESSBAR_FILE_NAME,
                                          SEND_EMAIL_WHEN_JOB_FINISHED_FROM_PIPELINE,
                                          CLEAN_OLD_JOBS_DIRECTORIES_FROM_PIPELINE)
from flask.SharedConsts import USER_FILE_NAME_ZIP, USER_FILE_NAME_TAR

PIPELINE_STEPS = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']


@dataclass
class Config:
    # Directories paths and file names - will be set by the pipeline
    run_dir: Path

    raw_data_path: Path
    data_path: Path
    genomes_names_path: Path

    steps_results_dir: Path
    tmp_dir: Path
    done_files_dir: Path
    final_output_dir: Path
    final_output_dir_name: str

    progressbar_file_path: Path
    error_file_path: Path

    # Job submission
    queue_name: str
    account_name: str
    node_name: str
    max_parallel_jobs: int
    use_job_manager: bool

    # Email at job completion
    email: str
    job_name: str
    run_number: str
    send_email: bool

    # Pipeline parameters
    identity_cutoff: float
    coverage_cutoff: float
    e_value_cutoff: float
    sensitivity: float
    core_minimal_percentage: float
    inputs_fasta_type: str

    outgroup: str
    bootstrap: bool
    max_number_of_core_ogs_for_phylogeny: int

    add_orphan_genes_to_ogs: bool
    run_optimized_mmseqs: bool
    filter_out_plasmids: bool
    num_of_genomes_in_batch: int
    pseudo_genome_mode: str

    # Debugging parameters
    do_not_copy_outputs_to_final_results_dir: bool
    bypass_number_of_genomes_limit: bool
    step_to_complete: str
    qfo_benchmark: bool
    only_calc_ogs: bool
    always_run_full_orthogroups_inference: bool
    use_parquet: bool
    verbose: bool

    # Clean at end
    clean_intermediate_outputs: bool
    clean_old_job_directories: bool

    # cpus, memory, time
    job_default_memory: str
    mmseqs_big_dataset_cpus: int
    mmseqs_big_dataset_memory: str
    mmseqs_time_limit: str
    phylogeny_cpus: int
    phylogeny_memory: str
    kegg_cpus: int
    kegg_memory: str
    codon_bias_cpus: int
    ani_cpus: int
    ani_memory: str
    orthoxml_memory: str
    phylogeny_time_limit: str
    infer_orthogroups_time_limit: str
    merge_sub_orthogroups_memory: str

    def to_csv(self, path: Path):
        config_df = pd.DataFrame(list(asdict(self).items()), columns=['key', 'value'])
        config_df.to_csv(path, index=False)

    @classmethod
    def from_csv(cls, path: Path):
        config_df = pd.read_csv(path, na_filter=False)
        data_dict = dict(zip(config_df["key"], config_df["value"]))

        typed_data = {}
        for key, value in data_dict.items():
            if key in cls.__annotations__:
                expected_type = cls.__annotations__[key]

                if expected_type == bool:
                    typed_data[key] = str_to_bool(value)
                else:
                    typed_data[key] = expected_type(value)
            else:
                raise ValueError(f'key {key} not in {cls} annotations')

        return cls(**typed_data)


def get_configuration():
    parser = argparse.ArgumentParser()
    parser.add_argument('--args_json_path', help='path to a json file that contains values for arguments which will '
                                                 'override the default values. Optional.')
    parser.add_argument('--run_dir', help='path to a directory where the pipeline will be run. Should contain a zip of'
                                          'the genomes. Mutually exclusive with --contigs_dir.')
    parser.add_argument('--contigs_dir',
                        help='path to a folder with the genomic sequences. This folder may be zipped, as well the files'
                             ' in it. Mutually exclusive with --run_dir.')
    parser.add_argument('--output_dir', help='relative path of directory where the output files will be written to',
                        default='outputs')
    parser.add_argument('--email',
                        help='A notification will be sent once the pipeline is done to this email. Optional.')
    parser.add_argument('--job_name', help='Optional job name.')
    parser.add_argument('--identity_cutoff', type=float, default=consts.DEFAULT_IDENTITY_CUTOFF,
                        help='minimum required percent of identity level (lower values will be filtered out)')
    parser.add_argument('--coverage_cutoff', type=float, default=consts.DEFAULT_COVERAGE_CUTOFF,
                        help='minimum required coverage percent of homology region between genes (lower values will be filtered out)')
    parser.add_argument('--e_value_cutoff', type=float, default=consts.DEFAULT_EVALUE_CUTOFF,
                        help='maxmimum permitted e-value (0 <= e_value_cutoff <= 1; higher values will be filtered'
                             ' out).')
    parser.add_argument('--sensitivity', type=float, default=consts.DEFAULT_MMSEQS_SENSITIVITY)
    parser.add_argument('--core_minimal_percentage', type=float, default=consts.DEFAULT_CORE_MINIMAL_PERCENTAGE,
                        help='the minimum required percent of gene members that is needed to be considered a core gene.'
                             ' For example: (1) 100 means that for a gene to be considered core, all strains should '
                             'have a member in the group.\n(2) 50 means that for a gene to be considered core, at least'
                             ' half of the strains should have a member in the group.\n(3) 0 means that every gene '
                             'should be considered as a core gene.')
    parser.add_argument('--bootstrap', type=str_to_bool, default=False,
                        help='whether or not to apply bootstrap procedure over the reconstructed species tree.')
    parser.add_argument('--outgroup',
                        help='The species name used to to root the phylogenetic tree, or empty string to leave unrooted.')
    parser.add_argument('--filter_out_plasmids', type=str_to_bool, default=False,
                        help='whether or not to filter out plasmids from the input files')
    parser.add_argument('--inputs_fasta_type', choices=['genomes', 'orfs'], default='genomes',
                        help='whether the input files are fastas of orfs and therefore Prodigal is skipped, '
                             'or genomes assemblies and then the first step is ORFs extraction using Prodigal')
    parser.add_argument('--add_orphan_genes_to_ogs', type=str_to_bool, default=False,
                        help='whether orphan genes should be considered as OGs')
    parser.add_argument('--qfo_benchmark', type=str_to_bool, default=False,
                        help='whether the input files are annotated genomes in the QfO benchmark format')
    # choices=['pupkoweb', 'pupkowebr', 'pupkolab', 'pupkolabr', 'pupkotmp', 'pupkotmpr', 'itaym', 'lilach',
    # 'bioseq', 'bental', 'oren.q', 'bioseq20.q'])
    parser.add_argument('-q', '--queue_name', help='The queue to which the job(s) will be submitted to',
                        default=consts.SLURM_PARTITION)
    parser.add_argument('--account_name', help='The slurm account to submit jobs to',
                        default=consts.SLURM_ACCOUNT)
    parser.add_argument('--node_name', help='The node name to submit jobs to')
    parser.add_argument('--step_to_complete', help='The final step to execute',
                        choices=[*PIPELINE_STEPS, ''])
    parser.add_argument('--only_calc_ogs', type=str_to_bool, default=False,
                        help='Do only the necessary steps to calculate OGs')
    parser.add_argument('--bypass_number_of_genomes_limit', type=str_to_bool, default=False,
                        help='Bypass the limit on number of genomes')
    parser.add_argument('--run_optimized_mmseqs', type=str_to_bool, default=False,
                        help='Optimize the mmseqs run')
    parser.add_argument('--use_parquet', type=str_to_bool, default=False,
                        help='When True, use parquet files when possible instead of csv')
    parser.add_argument('--do_not_copy_outputs_to_final_results_dir', type=str_to_bool, default=False, )
    parser.add_argument('--clean_intermediate_outputs', type=str_to_bool, default=False, )
    parser.add_argument('--num_of_genomes_in_batch', type=int, default=50)
    parser.add_argument('--pseudo_genome_mode', type=str, choices=['first_gene', 'consensus_gene'],
                        default='first_gene')
    parser.add_argument('--always_run_full_orthogroups_inference', type=str_to_bool, default=False, )
    parser.add_argument('--max_parallel_jobs', help='', type=int)
    parser.add_argument('--use_job_manager', type=str_to_bool, default=consts.USE_JOB_MANAGER)
    parser.add_argument('-v', '--verbose', type=str_to_bool, default=False,
                        help='Increase output verbosity')

    args = parser.parse_args()

    # Override arguments with args_json_path content
    if args.args_json_path:
        with open(args.args_json_path, 'r') as args_json_file:
            args_json = json.load(args_json_file)
            args.__dict__.update(args_json)

    # Validate arguments
    validate_arguments(args)

    (logger, times_logger, run_dir, error_file_path, progressbar_file_path, run_number, output_dir, tmp_dir,
     done_files_dir, steps_results_dir, data_path, final_output_dir_name, final_output_dir, genomes_names_path,
     raw_data_path) = \
        prepare_pipeline_framework(args)

    # Create Config object
    config = Config(run_dir=run_dir, raw_data_path=raw_data_path, data_path=data_path,
                    genomes_names_path=genomes_names_path,
                    steps_results_dir=steps_results_dir, tmp_dir=tmp_dir, done_files_dir=done_files_dir,
                    final_output_dir=final_output_dir, final_output_dir_name=final_output_dir_name,
                    error_file_path=error_file_path, progressbar_file_path=progressbar_file_path,

                    queue_name=args.queue_name, account_name=args.account_name, node_name=args.node_name,
                    max_parallel_jobs=args.max_parallel_jobs, use_job_manager=args.use_job_manager,

                    email=args.email, job_name=args.job_name, run_number=run_number,
                    send_email=consts.ENV == 'lsweb' and SEND_EMAIL_WHEN_JOB_FINISHED_FROM_PIPELINE and
                               not args.step_to_complete and not args.only_calc_ogs and
                               not args.do_not_copy_outputs_to_final_results_dir,

                    identity_cutoff=args.identity_cutoff, coverage_cutoff=args.coverage_cutoff,
                    e_value_cutoff=args.e_value_cutoff, sensitivity=args.sensitivity,
                    core_minimal_percentage=args.core_minimal_percentage, inputs_fasta_type=args.inputs_fasta_type,
                    outgroup=args.outgroup, bootstrap=args.bootstrap,
                    max_number_of_core_ogs_for_phylogeny=consts.MAX_NUMBER_OF_CORE_OGS_FOR_PHYLOGENY,

                    add_orphan_genes_to_ogs=args.add_orphan_genes_to_ogs,
                    filter_out_plasmids=args.filter_out_plasmids,
                    num_of_genomes_in_batch=args.num_of_genomes_in_batch, pseudo_genome_mode=args.pseudo_genome_mode,
                    run_optimized_mmseqs=args.run_optimized_mmseqs,

                    do_not_copy_outputs_to_final_results_dir=args.do_not_copy_outputs_to_final_results_dir,
                    bypass_number_of_genomes_limit=args.bypass_number_of_genomes_limit,
                    step_to_complete=args.step_to_complete,
                    qfo_benchmark=args.qfo_benchmark, only_calc_ogs=args.only_calc_ogs,
                    use_parquet=args.use_parquet,
                    always_run_full_orthogroups_inference=args.always_run_full_orthogroups_inference,
                    verbose=args.verbose,

                    clean_intermediate_outputs=args.clean_intermediate_outputs,
                    clean_old_job_directories=consts.ENV == 'lsweb' and CLEAN_OLD_JOBS_DIRECTORIES_FROM_PIPELINE,

                    job_default_memory=consts.JOB_DEFAULT_MEMORY_GB,
                    mmseqs_big_dataset_cpus=min(consts.MMSEQS_BIG_DATASET_CPUS, args.max_parallel_jobs),
                    mmseqs_big_dataset_memory=consts.MMSEQS_BIG_DATASET_MEMORY_GB,
                    mmseqs_time_limit=consts.MMSEQS_TIME_LIMIT,
                    phylogeny_cpus=min(consts.PHYLOGENY_CPUS, args.max_parallel_jobs),
                    phylogeny_memory=consts.PHYLOGENY_MEMORY_GB,
                    phylogeny_time_limit=consts.PHYLOGENY_TIME_LIMIT,
                    kegg_cpus=min(consts.KEGG_CPUS, args.max_parallel_jobs),
                    kegg_memory=consts.KEGG_MEMORY_GB,
                    codon_bias_cpus=min(consts.CODON_BIAS_CPUS, args.max_parallel_jobs),
                    ani_cpus=min(consts.ANI_CPUS, args.max_parallel_jobs),
                    ani_memory=consts.ANI_MEMORY_GB,
                    orthoxml_memory=consts.ORTHOXML_MEMORY_GB,
                    infer_orthogroups_time_limit=consts.INFER_ORTHOGROUPS_TIME_LIMIT,
                    merge_sub_orthogroups_memory=consts.MERGE_SUB_ORTHOGROUPS_MEMORY_GB
                    )

    config.to_csv(output_dir / 'config.csv')

    return logger, times_logger, config


def validate_arguments(args):
    if not args.max_parallel_jobs:
        if args.use_job_manager:
            args.max_parallel_jobs = consts.MAX_PARALLEL_JOBS
        else:
            args.max_parallel_jobs = os.cpu_count()

    if args.outgroup == "No outgroup":
        args.outgroup = ''

    if args.identity_cutoff < 0 or args.identity_cutoff > 100:
        raise ValueError(f'identity_cutoff argument {args.identity_cutoff} has invalid value')

    if args.coverage_cutoff < 0 or args.coverage_cutoff > 100:
        raise ValueError(f'coverage_cutoff argument {args.coverage_cutoff} has invalid value')

    if args.e_value_cutoff < 0 or args.e_value_cutoff > 1:
        raise ValueError(f'e_value_cutoff argument {args.e_value_cutoff} has invalid value')

    if args.core_minimal_percentage < 0 or args.core_minimal_percentage > 100:
        raise ValueError(f'core_minimal_percentage argument {args.core_minimal_percentage} has invalid value')


def prepare_pipeline_framework(args):
    if (args.run_dir and args.contigs_dir) or (not args.run_dir and not args.contigs_dir):
        raise ValueError('Either run_dir or contigs_dir should be provided, but not both.')

    if args.run_dir:
        run_dir = Path(args.run_dir)
        if (run_dir / USER_FILE_NAME_ZIP).exists():
            raw_data_path = run_dir / USER_FILE_NAME_ZIP
        elif (run_dir / USER_FILE_NAME_TAR).exists():
            raw_data_path = run_dir / USER_FILE_NAME_TAR
        else:
            raise ValueError(f'No genomes zip or tar file found in {run_dir}')
    else:  # args.contigs_dir was provided
        raw_data_path = Path(args.contigs_dir)
        if raw_data_path.exists():
            run_dir = raw_data_path.parent
        else:
            raise ValueError(f'contigs_dir argument {raw_data_path} does not exist!')

    output_dir = run_dir / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    level = logging.DEBUG if args.verbose else logging.INFO
    formatter = logging.Formatter(consts.LOG_MESSAGE_FORMAT)

    logger = logging.getLogger('main')
    main_file_handler = logging.FileHandler(output_dir / 'main_log.txt', mode='a')
    main_file_handler.setFormatter(formatter)
    logger.addHandler(main_file_handler)
    logger.setLevel(level)

    times_logger = logging.getLogger('times')
    times_file_handler = logging.FileHandler(output_dir / 'times_log.txt', mode='a')
    times_file_handler.setFormatter(formatter)
    times_logger.addHandler(times_file_handler)
    times_logger.setLevel(level)

    logger.info(args)
    logger.info(f'run_dir is: {run_dir}')
    logger.info(f'Created output_dir in: {output_dir}')

    final_output_dir_name = f'{WEBSERVER_NAME}_{args.output_dir}'
    final_output_dir = run_dir / final_output_dir_name
    logger.info(f'Creating final_output_dir is: {final_output_dir}')
    final_output_dir.mkdir(parents=True, exist_ok=True)

    error_file_path = final_output_dir / ERROR_FILE_NAME
    progressbar_file_path = run_dir / PROGRESSBAR_FILE_NAME

    run_number = run_dir.name
    logger.info(f'run_number is {run_number}')

    tmp_dir = output_dir / 'tmp'
    logger.info(f'Creating tmp_dir in: {tmp_dir}')
    tmp_dir.mkdir(parents=True, exist_ok=True)

    done_files_dir = output_dir / 'done'
    logger.info(f'Creating done_files_dir in: {done_files_dir}')
    done_files_dir.mkdir(parents=True, exist_ok=True)

    steps_results_dir = output_dir / 'steps_results'
    logger.info(f'Creating results_dir in: {steps_results_dir}')
    steps_results_dir.mkdir(parents=True, exist_ok=True)

    data_path = output_dir / 'inputs'
    data_path.mkdir(parents=True, exist_ok=True)

    genomes_names_path = output_dir / 'genomes_names.txt'

    return (logger, times_logger, run_dir, error_file_path, progressbar_file_path, run_number, output_dir, tmp_dir,
            done_files_dir, steps_results_dir, data_path, final_output_dir_name, final_output_dir,
            genomes_names_path, raw_data_path)
