# try:
import argparse
import subprocess
import sys
import filecmp
from time import ctime

sys.path.insert(0,"/groups/pupko/orenavr2/microbializer/") #ADD auxiliaries path
sys.path.insert(0,"/groups/pupko/orenavr2/microbializer/auxiliaries/") #ADD file_writer
from auxiliaries import file_writer #ADD file_writer
from auxiliaries.pipeline_auxiliaries import *

start = time()
logger = logging.getLogger('main')  # use logger instead of printing
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('contigs_dir', help='file containing commands for queue job ends with \t and a job name',
                    type=lambda file_path: str(file_path) if os.path.exists(file_path) else parser.error(f'{file_path} does not exist!'))
parser.add_argument('output_dir', help='directory where the log files will be written to')
parser.add_argument('email', help='A notification will be sent once the pipeline is done',
                    default='orenavram@gmail.com')
parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                    choices=['pupko', 'itaym', 'lilach', 'bioseq', 'bental'], default='pupko')
parser.add_argument('--dummy_delimiter',
                    help='The queue does not "like" very long commands. A dummy delimiter is used to break each row into different commands of a single job',
                    default='!@#')
parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
parser.add_argument('--src_dir', help='source code directory', type=lambda s: s.rstrip('/'), default='/groups/pupko/orenavr2/microbializer/pipeline/')

args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.DEBUG)

logger.fatal(args)

tmp_dir = os.path.join(args.output_dir, 'tmp_dir')
create_dir(tmp_dir)

done_files_dir = os.path.join(args.output_dir, 'done')
create_dir(done_files_dir)

data_path = os.listdir(args.contigs_dir)[0]
logger.info(data_path)
#if not os.path.isdir(data_path):
#    os.system(f'tar -xvf {data_path}')
#    data_path = os.path.join(data_path, os.listdir(data_path)[0])

# 1.	extract_orfs.py
# Input: (1) an input path for a fasta file with contigs/full genome (2) an output file path (with a suffix as follows: i_genes.fasta. especially relevant for the wrapper).
# Output: a fasta file where under each header, there’s a single gene.
# Can be parallelized on cluster
logger.warning(f'1: {"_"*80}')
dir_name = 'ORFs'
script_path = os.path.join(args.src_dir, 'extract_orfs.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '1_extract_orfs.txt')
if not os.path.exists(done_file_path):
    logger.warning('Extracting ORFs...')
    for fasta_file in os.listdir(args.contigs_dir):
        fasta_file_prefix = os.path.splitext(fasta_file)[0]
        output_file_name = f'{fasta_file_prefix}.{dir_name}'
        output_coord_name = f'{fasta_file_prefix}.gene_coordinates'
        params = [os.path.join(args.contigs_dir, fasta_file), os.path.join(pipeline_step_output_dir, output_file_name),
                  os.path.join(pipeline_step_output_dir, output_coord_name)] #Shir - path to translated sequences file
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules=['prodigal/prodigal-2.6.3'])
        num_of_expected_results += 1
    wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 2.  create_blast_DB.py
# Input: path to gene file to create DB from
# Output: Blast DB of the input file
# Can be parallelized on cluster
logger.warning(f'2: {"_"*80}')
ORFs_dir = previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'blast_db'
script_path = os.path.join(args.src_dir, 'create_blast_DB.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '2_create_blast_DB.txt')
if not os.path.exists(done_file_path):
    logger.warning('Creating Blast DB...')
    for fasta_file in os.listdir(previous_pipeline_step_output_dir):
        file_path = os.path.join(previous_pipeline_step_output_dir, fasta_file)
        fasta_file_prefix = fasta_file.split('.')[0]
        output_file_name = f'{fasta_file_prefix}.{dir_name}'
        params = [file_path, os.path.join(pipeline_step_output_dir,output_file_name)]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=fasta_file_prefix,queue_name=args.queue_name, required_modules=['blast/blast-2.2.30'])
        num_of_expected_results += 1
    wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 2.5.	blast_all_vs_all.py
# Input: (1) 2 input paths for 2 (different) genes files, g1 and g2 (2) an output file path (with a suffix as follows: i_vs_j_blast.tsv. especially relevant for the wrapper).
# Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
# Precisely, for each gene x in g1, blast x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
# Can be parallelized on cluster
logger.warning(f'2.5: {"_"*80}')
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'blast'
script_path = os.path.join(args.src_dir, 'blast_all_vs_all.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '2.5_blast_all_vs_all.txt')
if not os.path.exists(done_file_path):
    logger.warning('Blasting...')
    blasted_pairs = set()
    for fasta_file_name in os.listdir(ORFs_dir):
        fasta_file_path = os.path.join(ORFs_dir, fasta_file_name)
        strain1_name = fasta_file_name.split('.')[0]
        for db_file in os.listdir(previous_pipeline_step_output_dir):
            strain2_name = db_file.split('.')[0]
            pair = (strain1_name, strain2_name)
            if pair not in blasted_pairs and strain1_name != strain2_name:
                blasted_pairs.add(pair)
                logger.debug(f'{"#"*80}\nBlasting pair: {strain1_name}, {strain2_name}')
                logger.debug(blasted_pairs)
                db_prefix = os.path.splitext(os.path.join(previous_pipeline_step_output_dir, db_file))[0]
                output_file_name = f'{strain1_name}_vs_{strain2_name}.blast'
                params = [fasta_file_path, db_prefix, os.path.join(pipeline_step_output_dir, output_file_name)]
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules=['blast/blast-2.2.30'])
                num_of_expected_results += 1
    wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 3.	filter_blast.py
# Input: (1) a path for a i_vs_j_blast.tsv file (2) an output path (with a suffix as follows: i_vs_j_filtered.tsv. especially relevant for the wrapper).
# Output: the same format of the input file containing only pairs that passed the filtration. For each row in the input file (pair of genes), apply the following filters:
# 1. at least X% similarity
# 2. at least X% of the length
# 3.# write each pair to the output file if it passed all the above filters.
# Can be parallelized on cluster
logger.warning(f'3: {"_"*80}')
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'blast_filtered'
script_path = os.path.join(args.src_dir, 'filter_blast_results.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '3_filter_blast_results.txt')
if not os.path.exists(done_file_path):
    logger.warning('Filtering...')
    for blast_results_file in os.listdir(previous_pipeline_step_output_dir):
        fasta_file_prefix = os.path.splitext(blast_results_file)[0]
        output_file_name = f'{fasta_file_prefix}.{dir_name}'
        params = [os.path.join(previous_pipeline_step_output_dir, blast_results_file), os.path.join(pipeline_step_output_dir, output_file_name)]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
        num_of_expected_results += 1
    wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 4.	find_reciprocal_hits.py
# Input: (1) a path for a i_vs_j_filtered.tsv file and a path for a j_vs_i_filtered.tsv file (2) an output path (with a suffix as follows: i_vs_j_reciprocal_hits.tsv
# Output: a tab delimited file containing only reciprocal best-hit pairs and their bit score from blast, i.e., if x’s best hit was y in the first file with bit score z, x\ty\tz will appear in the output file only if y’s best hit in the other file will be x.
# Can be parallelized on cluster
logger.warning(f'4: {"_"*80}')
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'reciprocal_hits'
script_path = os.path.join(args.src_dir, 'reciprocal_hits.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '4_reciprocal_hits.txt')
if not os.path.exists(done_file_path):
    logger.warning('Extracting reciprocal hits...')
    for blast_filtered_results_file in os.listdir(previous_pipeline_step_output_dir):
        fasta_file_prefix = os.path.splitext(blast_filtered_results_file)[0]
        id_1 = fasta_file_prefix.split("_vs_")[0]
        id_2 = fasta_file_prefix.split("_vs_")[1]
        if (id_1<id_2):
            reciprocal_file_prefix = "_vs_".join(fasta_file_prefix.split("_vs_")[::-1])
            reciprocal_file_name = reciprocal_file_prefix + "." + blast_filtered_results_file.split('.')[1]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            params = [os.path.join(previous_pipeline_step_output_dir, blast_filtered_results_file), os.path.join(previous_pipeline_step_output_dir, reciprocal_file_name), pipeline_step_output_dir, output_file_name]
            print("\n\n\n",script_path,params[0], params[1] , params[2] ,"\n\n\n\n")
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
            num_of_expected_results += 1
    wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 4.5. concatenate_reciprocal_hits
# Input: path to folder with all reciprocal hits files
# Output: concatenated file of all reciprocal hits files
# CANNOT be parallelized on cluster
logger.warning(f'4.5: {"_"*80}')
done_file_path = os.path.join(done_files_dir, '4.5_concatenate_reciprocal_hits.txt')
if not os.path.exists(done_file_path):
    logger.warning('Concatenating reciprocal hits...')
    cmd = f'cat {pipeline_step_output_dir}/*.{dir_name} > {args.output_dir}/all_recip_hits.csv'
    subprocess.call(cmd, shell = True)
    # No need to wait...
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 5.	construct_putative_orthologs_table.py
# Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
# Output: updates the table with the info from the reciprocal hit file.
# CANNOT be parallelized on cluster
logger.warning(f'5: {"_"*80}')
previous_pipeline_step_output_dir = pipeline_step_output_dir
putative_orthologs_table_path = os.path.join(args.output_dir, 'putative_orthologs_table.tsv')
script_path = os.path.join(args.src_dir, 'construct_putative_orthologs_table.py')
done_file_path = os.path.join(done_files_dir, '5_construct_putative_orthologs_table.txt')
if not os.path.exists(done_file_path):
    logger.warning('Constructing putative orthologs table...')
    for reciprocal_hits_file in os.listdir(previous_pipeline_step_output_dir):
        subprocess.check_output(['python', '-u', script_path, os.path.join(previous_pipeline_step_output_dir, reciprocal_hits_file), putative_orthologs_table_path])
    # No need to wait...
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 6.	split_putative_orthologs_table
# Input: (1) a path for a putative orthologs table (each row in this table is a putative orthologous set) (2) an output_path to a directory.
# Output: each line is written to a separate file in the output directory.
# Can be parallelized on cluster (better as a subprocess)
# logger.warning('Splitting putative orthologs table...')
logger.warning(f'6: {"_"*80}')
dir_name = 'splitted_putative_orthologs_table'
script_path = os.path.join(args.src_dir, '..', 'auxiliaries', 'file_writer.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '6_split_putative_orthologs_table.txt')
if not os.path.exists(done_file_path):
    logger.warning('Splitting putative orthologs table...')
    with open(putative_orthologs_table_path) as f:
        f.readline()  # skip header
        lines = ''
        i = 1
        for line in f:
            out_file = os.path.join(pipeline_step_output_dir, (line.split(',')[0]).split('|')[1] + '.' + dir_name)
            lines+=line
            #subprocess.call(['python', script_path, out_file,'--content' ,line])
            if (i % 100 == 0):
                file_writer.write_to_file(out_file, lines)
                lines = ''
            #num_of_expected_results += 1
            i += 1
    if lines:
        file_writer.write_to_file(out_file, lines)
    # wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 7.	verify_cluster.py
# Input: (1) a file path (containing a row from the putative orthologs table) (2) an output path to a directory.
# Output: copies this file to the output folder if the group is well-clustered.
# Can be parallelized on cluster
logger.warning(f'7: {"_"*80}')
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'splitted_final_orthologs_table'
script_path = os.path.join(args.src_dir, 'verify_cluster.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '7_verify_cluster.txt')
if not os.path.exists(done_file_path):
    logger.warning('Verifying clusters...')
    for putative_orthologs_set in os.listdir(previous_pipeline_step_output_dir):
        putative_orthologs_set_prefix = os.path.splitext(putative_orthologs_set)[0]
        output_file_name = f'{putative_orthologs_set_prefix}.{dir_name}'
        params = [os.path.join(args.output_dir, 'all_recip_hits.csv'),
                  os.path.join(previous_pipeline_step_output_dir, putative_orthologs_set), pipeline_step_tmp_dir,
                  os.path.join(pipeline_step_output_dir, output_file_name)]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules=['MCL-edge/mcl-14-137'])
        num_of_expected_results += 1
    wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, time_to_wait=5)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 8.	Join_final_orthologs_table.py
# Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
# Output: aggregates all the well-clustered OGs to the final table.
logger.warning(f'8: {"_"*80}')
previous_pipeline_step_output_dir = pipeline_step_output_dir
final_orthologs_table_path = os.path.join(args.output_dir, 'final_orthologs_table.tsv')
script_path = os.path.join(args.src_dir, 'construct_final_table.py')
done_file_path = os.path.join(done_files_dir, '8_construct_final_table.txt')
if not os.path.exists(done_file_path):
    logger.warning('Joining clusters to final orthologs table...')
    with open(putative_orthologs_table_path) as f:
        orthologs_table_header = f.readline().rstrip()
    with open(final_orthologs_table_path, 'w') as f:
        f.write(orthologs_table_header + '\n')
    for final_orthologs_set_file in os.listdir(previous_pipeline_step_output_dir):
        cmd = f'cat {os.path.join(previous_pipeline_step_output_dir, final_orthologs_set_file)} >> {final_orthologs_table_path}'
        subprocess.check_output(cmd, shell=True)
    # No need to wait...
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')


# 9.	extract_orthologs_sequences.py
# Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
# Output: write the sequences of the orthologs group to the output file.
logger.warning(f'9: {"_"*80}')
dir_name = 'orthologs_sets_sequences'
script_path = os.path.join(args.src_dir, 'extract_orthologs_sequences.py')
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
done_file_path = os.path.join(done_files_dir, '9_extract_orthologs_sequences.txt')
if not os.path.exists(done_file_path):
    logger.warning('Extracting orthologs set sequences according to final orthologs table...')
    for final_orthologs_set_file in os.listdir(previous_pipeline_step_output_dir):
        fasta_file_prefix = os.path.splitext(final_orthologs_set_file)[0]
        output_file_name = f'{fasta_file_prefix}.{dir_name}'
        params = [os.path.join(previous_pipeline_step_output_dir, final_orthologs_set_file), os.path.join(previous_pipeline_step_output_dir, 'ORFs'), pipeline_step_output_dir]
        print(f'params for orthologs_sets_sequences: {params}')
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
        num_of_expected_results += 1
    wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results)
    file_writer.write_to_file(done_file_path)
else:
    logger.warning(f'done file {done_file_path} already exists. Skipping step...')

status = 'is done'

# except:
#     import logging
#     logger = logging.getLogger('main')  # use logger instead of printing
#
#     status = 'was failed'
#     exc_type, exc_obj, exc_tb = sys.exc_info()
#     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#     logger.error('\n\n' + '$' * 60 + '\n\n')
#     logger.error(f'{ctime()}: microbializer failed :(\n\n')
#     logger.error(f'{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\n')
#     logger.error('$' * 60)

end = time()
msg = f'microbializer pipeline {status}. Results can be found at {args.output_dir}. Took {measure_time(int(end-start))}'
logger.warning(msg)
send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>', args.email, subject=f'Microbialzer {status}.', content=msg)

