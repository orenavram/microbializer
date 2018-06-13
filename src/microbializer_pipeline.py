import os
import logging
import argparse
import sys
from time import time
from pipeline_auxiliaries import *

start = time()
logger = logging.getLogger('main')  # use logger instead of printing

parser = argparse.ArgumentParser()
parser.parse_args()

#TODO: get arguments from the user

contigs_dir = ''
output_dir = ''
raw_file_suffix = '.tsv'
src_dir = ''
new_line_delimiter = '!@#'
tmp_dir = os.path.join(output_dir, 'tmp_dir')
create_dir(tmp_dir)
# q_submitter_script_path = '/bioseq/bioSequence_scripts_and_constants/q_submitter.py'
# done_files_script_path = os.path.join(src_dir, 'file_writer.py')

queue_name = 'pupko'
# 1.	extract_orfs.py
# Input: (1) an input path for a fasta file with contigs/full genome (2) an output file path (with a suffix as follows: i_genes.fasta. especially relevant for the wrapper).
# Output: a fasta file where under each header, there’s a single gene.
# Can be parallelized on cluster
dir_name = 'ORFs'
script_name = 'extract_orfs.py'
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(output_dir, tmp_dir, dir_name)
for fasta_file in os.listdir(contigs_dir):
    fasta_file_prefix = fasta_file.split('.')[0]
    output_file_name = fasta_file_prefix + '.' + dir_name
    params = [os.path.join(contigs_dir, fasta_file), os.path.join(pipeline_step_output_dir, output_file_name)]
    submit_pipeline_step(script_name, params, pipeline_step_tmp_dir, job_name=output_file_name)
    num_of_expected_results += 1

wait_for_results(pipeline_step_tmp_dir, script_name, num_of_expected_results)


# 2.	blast_all_vs_all.py
# Input: (1) 2 input paths for 2 (different) genes files, g1 and g2 (2) an output file path (with a suffix as follows: i_vs_j_blast.tsv. especially relevant for the wrapper).
# Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
# Precisely, for each gene x in g1, blast x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
# Can be parallelized on cluster
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'blast'
script_name = 'blast_all_vs_all.py'
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(output_dir, tmp_dir, dir_name)
for fasta_file1 in os.listdir(previous_pipeline_step_output_dir):
    file1_path = os.path.join(previous_pipeline_step_output_dir, fasta_file1)
    fasta_file1_prefix = fasta_file1.split('.')[0]
    output_file_name = fasta_file1_prefix + '.blast'
    for fasta_file2 in os.listdir(previous_pipeline_step_output_dir):
        file2_path = os.path.join(previous_pipeline_step_output_dir, fasta_file2)
        file2_name_prefix = fasta_file2.split('.')[0]
        output_file_name = '_'.join([fasta_file1_prefix, 'vs', file2_name_prefix + '.blast'])
        params = [file1_path, file2_path, output_file_name]
        submit_pipeline_step(script_name, params, pipeline_step_tmp_dir, job_name=output_file_name)
        num_of_expected_results += 1

wait_for_results(pipeline_step_tmp_dir, script_name, num_of_expected_results)


# 3.	filter_blast.py
# Input: (1) a path for a i_vs_j_blast.tsv file (2) an output path (with a suffix as follows: i_vs_j_filtered.tsv. especially relevant for the wrapper).
# Output: the same format of the input file containing only pairs that passed the filtration. For each row in the input file (pair of genes), apply the following filters:
# 1. at least X% similarity
# 2. at least X% of the length
# 3.# write each pair to the output file if it passed all the above filters.
# Can be parallelized on cluster
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'blast_filtered'
script_name = 'filter_blast.py'
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(output_dir, tmp_dir, dir_name)
for blast_results_file in os.listdir(previous_pipeline_step_output_dir):
    fasta_file_prefix = blast_results_file.split('.')[0]
    output_file_name = fasta_file_prefix + '.' + dir_name
    params = [os.path.join(previous_pipeline_step_output_dir, blast_results_file), os.path.join(pipeline_step_output_dir, output_file_name)]
    submit_pipeline_step(script_name, params, pipeline_step_tmp_dir, job_name=output_file_name)
    num_of_expected_results += 1

wait_for_results(pipeline_step_tmp_dir, script_name, num_of_expected_results)


# 4.	find_reciprocal_hits.py
# Input: (1) a path for a i_vs_j_filtered.tsv file and a path for a j_vs_i_filtered.tsv file (2) an output path (with a suffix as follows: i_vs_j_reciprocal_hits.tsv
# Output: a tab delimited file containing only reciprocal best-hit pairs and their bit score from blast, i.e., if x’s best hit was y in the first file with bit score z, x\ty\tz will appear in the output file only if y’s best hit in the other file will be x.
# Can be parallelized on cluster
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'reciprocal_hits'
script_name = 'find_reciprocal_hits.py'
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(output_dir, tmp_dir, dir_name)
for blast_filtered_results_file in os.listdir(previous_pipeline_step_output_dir):
    fasta_file_prefix = blast_filtered_results_file.split('.')[0]
    output_file_name = fasta_file_prefix + '.' + dir_name
    params = [os.path.join(previous_pipeline_step_output_dir, blast_filtered_results_file), os.path.join(pipeline_step_output_dir, output_file_name)]
    submit_pipeline_step(script_name, params, pipeline_step_tmp_dir, job_name=output_file_name)
    num_of_expected_results += 1

wait_for_results(pipeline_step_tmp_dir, script_name, num_of_expected_results)


# 5.	construct_putative_orthologs_table.py
# Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
# Output: updates the table with the info from the reciprocal hit file.
# CANNOT be parallelized on cluster
previous_pipeline_step_output_dir = pipeline_step_output_dir
putative_orthologs_table_path = os.path.join(output_dir, 'putative_orthologs_table.tsv')
script_name = 'construct_putative_orthologs_table.py'
for reciprocal_hits_file in os.listdir(previous_pipeline_step_output_dir):
    os.system(f'{script_name} {os.path.join(previous_pipeline_step_output_dir, reciprocal_hits_file)} {putative_orthologs_table_path}')

# No need to wait...


# 6.	split_putative_orthologs_table.py
# Input: (1) a path for a putative orthologs table (each row in this table is a putative orthologous set) (2) an output_path to a directory.
# Output: each line is written to a separate file in the output directory.
# Can be parallelized on cluster (better as a subprocess)
dir_name = 'splitted_putative_orthologs_table'
script_name = 'split_putative_orthologs_table.py'
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(output_dir, tmp_dir, dir_name)
with open(putative_orthologs_table_path) as f:
    for line in f:
        out_file = os.path.join(pipeline_step_output_dir, line.split('\t')[0] + '.' + dir_name)
        subprocess.call(['python', 'file_writer.py', out_file, line])
    num_of_expected_results += 1

wait_for_results(pipeline_step_output_dir, 'splitting_putative_orthologs_table', num_of_expected_results)


# 7.	verify_cluster.py
# Input: (1) a file path (containing a row from the putative orthologs table) (2) an output path to a directory.
# Output: copies this file to the output folder if the group is well-clustered.
# Can be parallelized on cluster
previous_pipeline_step_output_dir = pipeline_step_output_dir
dir_name = 'splitted_final_orthologs_table'
script_name = 'verify_cluster.py'
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(output_dir, tmp_dir, dir_name)
for putative_orthologs_set in os.listdir(previous_pipeline_step_output_dir):
    putative_orthologs_set_prefix = putative_orthologs_set.split('.')[0]
    output_file_name = putative_orthologs_set_prefix + '.' + dir_name
    params = [os.path.join(previous_pipeline_step_output_dir, putative_orthologs_set), os.path.join(pipeline_step_output_dir, output_file_name)]
    submit_pipeline_step(script_name, params, pipeline_step_tmp_dir, job_name=output_file_name)
    num_of_expected_results += 1

wait_for_results(pipeline_step_tmp_dir, script_name, num_of_expected_results)


#TODO: from here and on
# 8.	Join_final_orthologs_table.py
# Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
# Output: aggregates all the well-clustered OGs to the final table.
previous_pipeline_step_output_dir = pipeline_step_output_dir
final_orthologs_table_path = os.path.join(output_dir, 'final_orthologs_table.tsv')
script_name = 'construct_final_table.py'
with open(putative_orthologs_table_path) as f:
    orthologs_table_header = f.readline().rstrip()
with open(final_orthologs_table_path, 'w') as f:
    f.write(orthologs_table_header + '\n')
    for final_orthologs_set_file_file in os.listdir(previous_pipeline_step_output_dir):
        os.system(f'cat {os.path.join(previous_pipeline_step_output_dir, final_orthologs_set_file_file)} >> {final_orthologs_table_path}')

# No need to wait...


# 9.	extract_sequences.py
# Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
# Output: write the sequences of the orthologs group to the output file.
dir_name = 'orthologs_sets_sequences'
script_name = 'extract_sequences.py'
num_of_expected_results = 0
pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(output_dir, tmp_dir, dir_name)
for final_orthologs_set_file_file in os.listdir(previous_pipeline_step_output_dir):
    fasta_file_prefix = final_orthologs_set_file_file.split('.')[0]
    output_file_name = fasta_file_prefix + '.' + dir_name
    params = [os.path.join(previous_pipeline_step_output_dir, final_orthologs_set_file_file), pipeline_step_output_dir]
    submit_pipeline_step(script_name, params, pipeline_step_tmp_dir, job_name=output_file_name)
    num_of_expected_results += 1

wait_for_results(pipeline_step_tmp_dir, script_name, num_of_expected_results)

end = time()
logger.info(f'Microbializer pipeline is done. Results can be found at {output_dir}. Took {measure_time(start, end)}')