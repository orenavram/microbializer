"""
script_name.py /Users/Oren/Dropbox/Projects/microbializer/data_for_test_cds/ /Users/Oren/Dropbox/Projects/microbializer/output_examples/mock_output/ orenavram@gmail.com -q pupko
"""


def notify_admin(meta_output_dir, meta_output_url, run_number, CONSTS):
    email = 'NO_EMAIL'
    user_email_path = os.path.join(meta_output_dir, 'user_email.txt')
    if os.path.exists(user_email_path):
        with open(user_email_path) as f:
            email = f.read().rstrip()
    error_log_path = 'NO_ERROR_LOG'
    tmp = [file for file in os.listdir(meta_output_dir) if file.endswith('.err')]
    if len(tmp) > 0:
        error_log_path = tmp[0]
    # Send me a notification email every time there's a failure
    send_email(smtp_server=CONSTS.SMTP_SERVER,
               sender=CONSTS.ADMIN_EMAIL,
               receiver=CONSTS.OWNER_EMAIL,
               subject=f'{CONSTS.WEBSERVER_NAME} job {run_number} by {email} has been failed: ',
               content=f"{email}\n\n{os.path.join(meta_output_url, 'output.html')}\n\n"
               f"{os.path.join(meta_output_url, 'cgi_debug.txt')}\n\n"
               f"{os.path.join(meta_output_url, error_log_path)}\n\n"
               f"{os.path.join(meta_output_dir, error_log_path)}")


def add_results_to_final_dir(source, final_output_dir):
    dest = os.path.join(final_output_dir, os.path.split(source)[1])
    logger.info(f'Moving {source} TO {dest}')
    # shutil.move(source, dest)
    try:
        shutil.copytree(source, dest)
    except FileExistsError:
        pass
    return dest


try:
    import argparse
    import sys
    import os
    import tarfile
    import shutil

    print(os.getcwd())
    print(f'sys.path is\n{sys.path}')

    if os.path.exists('/bioseq'):  # remote run
        remote_run = True
        sys.path.append('/bioseq/microbializer/auxiliaries')
        sys.path.append('/bioseq/microbializer/cgi')
        sys.path.append('/bioseq/bioSequence_scripts_and_constants/') #ADD file_writer
    else:
        # local run
        remote_run = False
        sys.path.append('../auxiliaries')
        sys.path.append('../cgi')

    print(f'sys.path is\n{sys.path}')

    import file_writer #ADD file_writer
    from email_sender import send_email
    from pipeline_auxiliaries import *
    from plots_generator import *

    import WEBSERVER_CONSTANTS as CONSTS

    from html_editor import edit_success_html, edit_failure_html

    start = time()

    parser = argparse.ArgumentParser()
    parser.add_argument('contigs_dir', help='path to a folder with the genomic sequences. This folder may be zipped, as well the files in it.',
                        type=lambda path: path.rstrip('/') if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('output_dir', help='directory where the output files will be written to',
                        type=lambda path: path.rstrip('/'))
    parser.add_argument('email', help='A notification will be sent once the pipeline is done',
                        default=CONSTS.OWNER_EMAIL)
    parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                        choices=['pupko', 'itaym', 'lilach', 'bioseq', 'bental', 'oren.q', 'bioseq20.q'], default='pupko')
    parser.add_argument('--dummy_delimiter',
                        help='The queue does not "like" very long commands. A dummy delimiter is used to break each row into different commands of a single job',
                        default='!@#')
    parser.add_argument('--src_dir', help='source code directory', type=lambda s: s.rstrip('/'), default='/bioseq/microbializer/pipeline')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    logger.info(args)

    create_dir(args.output_dir)

    run_number = args.output_dir.split('/')[-2]

    meta_output_dir = os.path.join(args.output_dir, '..')
    meta_output_url = os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number)

    output_html_path = os.path.join(meta_output_dir, 'output.html')
    logger.info(f'output_html_path is {output_html_path}')

    logger.info(f'run_number is {run_number}')

    output_url = os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number, 'output.html')
    logger.info(f'output_url is {output_url}')

    error_file_path = os.path.join(args.output_dir, 'error.txt')

    tmp_dir = os.path.join(args.output_dir, 'tmp_dir')
    create_dir(tmp_dir)

    done_files_dir = os.path.join(args.output_dir, 'done')
    create_dir(done_files_dir)

    data_path = args.contigs_dir
    logger.info(f'data_path is: {data_path}')


    # extract tar folder
    if not os.path.isdir(data_path):
        if tarfile.is_tarfile(data_path):
            with tarfile.open(data_path, 'r:gz') as f:
                f.extractall(path=meta_output_dir) # unzip tar folder to parent dir
            data_path = data_path.split('.tar')[0] # e.g., /groups/pupko/orenavr2/microbializer/example_data.tar.gz
            logger.info(f'Updated data_path is:\n{data_path}')
        else: #TODO: check if the else works!
            shutil.unpack_archive(data_path, extract_dir=meta_output_dir) # unzip tar folder to parent dir
            data_path = os.path.splitext(data_path)[0] # e.g., /groups/pupko/orenavr2/microbializer/example_data.tar.gz
            logger.info(f'Updated data_path is:\n{data_path}')

    assert os.path.isdir(data_path)

    # gunzip gz files in $data_path if any
    for file in os.listdir(data_path):
        #TODO: handle the case where the zipped item is a folder of a folder of a folder (etc...) of the genomes
        file_path = os.path.join(data_path, file)
        if file_path.endswith('gz'):
            subprocess.run(f'gunzip -f {file_path}', shell=True)


    # 1.	extract_orfs.py
    # Input: (1) an input path for a fasta file with contigs/full genome (2) an output file path (with a suffix as follows: i_genes.fasta. especially relevant for the wrapper).
    # Output: a fasta file where under each header, there’s a single gene.
    # Can be parallelized on cluster
    # Prodigal ignores newlines in the middle of a sequence, namely, >bac1\nAAA\nAA\nTTT >bac1\nAAAAATTT will be analyzed identically.
    step = '01'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_ORFs'
    script_path = os.path.join(args.src_dir, 'extract_orfs.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_extract_orfs.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting ORFs...')
        for fasta_file in os.listdir(data_path):
            fasta_file_prefix = os.path.splitext(fasta_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            output_coord_name = f'{fasta_file_prefix}.gene_coordinates'
            params = [os.path.join(data_path, fasta_file),
                      os.path.join(pipeline_step_output_dir, output_file_name),
                      os.path.join(pipeline_step_output_dir, output_coord_name)] #Shir - path to translated sequences file
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules_as_list=['prodigal/prodigal-2.6.3'])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 2.  create_blast_DB.py
    # Input: path to gene file to create DB from
    # Output: Blast DB of the input file
    # Can be parallelized on cluster
    step = '02'
    logger.info(f'Step {step}: {"_"*100}')
    ORFs_dir = previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_blast_db'
    script_path = os.path.join(args.src_dir, 'create_blast_DB.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_create_blast_DB.txt')
    if not os.path.exists(done_file_path):
        logger.info('Creating Blast DB...')
        for fasta_file in os.listdir(previous_pipeline_step_output_dir):
            file_path = os.path.join(previous_pipeline_step_output_dir, fasta_file)
            fasta_file_prefix = os.path.splitext(fasta_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            params = [file_path,
                      os.path.join(pipeline_step_output_dir,output_file_name)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=fasta_file_prefix,queue_name=args.queue_name, required_modules_as_list=['blast/blast-2.2.30'])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 3.	blast_all_vs_all.py
    # Input: (1) 2 input paths for 2 (different) genes files, g1 and g2 (2) an output file path (with a suffix as follows: i_vs_j_blast.tsv. especially relevant for the wrapper).
    # Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
    # Precisely, for each gene x in g1, blast x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
    # Can be parallelized on cluster
    step = '03'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_blast_analysis'
    script_path = os.path.join(args.src_dir, 'blast_all_vs_all.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_blast_all_vs_all.txt')
    if not os.path.exists(done_file_path):
        logger.info('Blasting...')
        blasted_pairs = set()
        for fasta_file_name in os.listdir(ORFs_dir):
            fasta_file_path = os.path.join(ORFs_dir, fasta_file_name)
            strain1_name = os.path.splitext(fasta_file_name)[0]
            for db_file in os.listdir(previous_pipeline_step_output_dir):
                # 2 extensions should be removed GCF_000006945.2_ASM694v2_cds_from_genomic.02_blast_db.nhr or GCF_000006945.2_ASM694v2_cds_from_genomic.02_blast_db.nin
                strain2_name = os.path.splitext(os.path.splitext(db_file)[0])[0]
                pair = (strain1_name, strain2_name)
                if pair not in blasted_pairs and strain1_name != strain2_name: #TODO: WHY THE HELL this condition is used? a potential bug!
                    blasted_pairs.add(pair)
                    logger.debug(f'{"#"*100}\nBlasting pair: {strain1_name}, {strain2_name}')
                    logger.debug(blasted_pairs)
                    db_prefix = os.path.splitext(os.path.join(previous_pipeline_step_output_dir, db_file))[0]
                    output_file_name = f'{strain1_name}_vs_{strain2_name}.{dir_name}'
                    params = [fasta_file_path,
                              db_prefix,
                              os.path.join(pipeline_step_output_dir, output_file_name)]
                    submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules_as_list=['blast/blast-2.2.30'])
                    num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 4.	filter_blast.py
    # Input: (1) a path for a i_vs_j_blast.tsv file (2) an output path (with a suffix as follows: i_vs_j_filtered.tsv. especially relevant for the wrapper).
    # Output: the same format of the input file containing only pairs that passed the filtration. For each row in the input file (pair of genes), apply the following filters:
    # 1. at least X% similarity
    # 2. at least X% of the length
    # 3.# write each pair to the output file if it passed all the above filters.
    # Can be parallelized on cluster
    step = '04'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_blast_filtered'
    script_path = os.path.join(args.src_dir, 'filter_blast_results.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_filter_blast_results.txt')
    if not os.path.exists(done_file_path):
        logger.info('Filtering...')
        for blast_results_file in os.listdir(previous_pipeline_step_output_dir):
            fasta_file_prefix = os.path.splitext(blast_results_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            params = [os.path.join(previous_pipeline_step_output_dir, blast_results_file),
                      os.path.join(pipeline_step_output_dir, output_file_name)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 5.	find_reciprocal_hits.py
    # Input: (1) a path for a i_vs_j.blast_filtered file and a path for a j_vs_i.blast_filtered file (2) an output path (with a suffix as follows: i_vs_j_reciprocal_hits.tsv
    # Output: a tab delimited file containing only reciprocal best-hit pairs and their bit score from blast, i.e., if x’s best hit was y in the first file with bit score z, x\ty\tz will appear in the output file only if y’s best hit in the other file will be x.
    # Can be parallelized on cluster
    step = '05'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_reciprocal_hits'
    script_path = os.path.join(args.src_dir, 'find_reciprocal_hits.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_find_reciprocal_hits.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting reciprocal hits...')
        for blast_filtered_results_file in os.listdir(previous_pipeline_step_output_dir):
            file_name_prefix = os.path.splitext(blast_filtered_results_file)[0]
            strain1, strain2 = file_name_prefix.split('_vs_')[:2]
            logger.debug(f'\nstrain1: {strain1}\nstrain2: {strain2}')
            if strain1<strain2:
                # avoid both strain1,strain2 and strain2,strain1
                reciprocal_file_name_prefix = "_vs_".join([strain2, strain1])
                reciprocal_file_name = reciprocal_file_name_prefix + os.path.splitext(blast_filtered_results_file)[1]
                output_file_name = f'{file_name_prefix}.{dir_name}'
                params = [os.path.join(previous_pipeline_step_output_dir, blast_filtered_results_file),
                          os.path.join(previous_pipeline_step_output_dir, reciprocal_file_name),
                          os.path.join(pipeline_step_output_dir, output_file_name)]
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
                num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 6. concatenate_reciprocal_hits
    # Input: path to folder with all reciprocal hits files
    # Output: concatenated file of all reciprocal hits files
    # CANNOT be parallelized on cluster
    step = '06'
    logger.info(f'Step {step}: {"_"*100}')
    all_reciprocal_hits_file = os.path.join(args.output_dir, 'concatenated_all_reciprocal_hits.txt')
    done_file_path = os.path.join(done_files_dir, f'{step}_concatenate_reciprocal_hits.txt')
    if not os.path.exists(done_file_path):
        logger.info('Concatenating reciprocal hits...')
        cmd = f'cat {pipeline_step_output_dir}/*.{dir_name} > {all_reciprocal_hits_file}'
        subprocess.run(cmd, shell = True)
        # No need to wait...
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 7.	construct_putative_orthologs_table.py
    # Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
    # Output: updates the table with the info from the reciprocal hit file.
    # CANNOT be parallelized on cluster
    step = '07'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_putative_table'
    script_path = os.path.join(args.src_dir, 'construct_putative_orthologs_table.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    putative_orthologs_table_path = os.path.join(pipeline_step_output_dir, 'putative_orthologs_table.txt')
    done_file_path = os.path.join(done_files_dir, f'{step}_construct_putative_orthologs_table.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing putative orthologs table...')
        job_name = os.path.split(script_path)[-1]
        params = [all_reciprocal_hits_file,
                  putative_orthologs_table_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name,
                             queue_name=args.queue_name)
        num_of_expected_results = 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # # 6.	split_putative_orthologs_table
    # # Input: (1) a path for a putative orthologs table (each row in this table is a putative orthologous group) (2) an output_path to a directory.
    # # Output: each line is written to a separate file in the output directory.
    # # Can be parallelized on cluster (better as a subprocess)
    # logger.info(f'Step 6: {"_"*100}')
    # dir_name = 'splitted_putative_orthologs_table'
    # # script_path = os.path.join(args.src_dir, '..', 'auxiliaries', 'file_writer.py')
    # # num_of_expected_results = 0
    # pipeline_step_output_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)[0]  # pipeline_step_tmp_dir is not needed
    # done_file_path = os.path.join(done_files_dir, '6_split_putative_orthologs_table.txt')
    # if not os.path.exists(done_file_path):
    #     logger.info('Splitting putative orthologs table...')
    #     with open(putative_orthologs_table_path) as f:
    #         header = f.readline()
    #         lines = ''
    #         i = 0
    #         for line in f:
    #             lines += line
    #             i += 1
    #             if i % 100 == 0:
    #                 out_file = os.path.join(pipeline_step_output_dir, (line.split(',')[0]) + '.' + dir_name)
    #                 # subprocess.call(['python', script_path, out_file,'--content' ,line])
    #                 file_writer.write_to_file(out_file, header+lines)
    #                 lines = ''
    #         if lines: # write last rows (maybe we didn't reach 100...)
    #             out_file = os.path.join(pipeline_step_output_dir, (line.split(',')[0]) + '.' + dir_name)
    #             file_writer.write_to_file(out_file, header+lines)
    #         # wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
    #
    #     file_writer.write_to_file(done_file_path)
    # else:
    #     logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 8   prepare_files_for_mcl.py
    # Input: (1) a path for a concatenated all reciprocal hits file (2) a path for a putative orthologs file (3) a path for an output folder
    # Output: an input file for MCL for each putative orthologs group
    # CANNOT be parallelized on cluster (if running on the concatenated file)
    step = '08'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_mcl_input_files'
    script_path = os.path.join(args.src_dir, 'prepare_files_for_mcl.py')
    num_of_expected_results = 1 # a single job that prepares all the files
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_prepare_files_for_mcl.txt')
    if not os.path.exists(done_file_path):
        logger.info('Preparing files for MCL...')
        job_name = os.path.split(script_path)[-1]
        params = [all_reciprocal_hits_file,
                  putative_orthologs_table_path,
                  pipeline_step_output_dir]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name, queue_name=args.queue_name)
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 9.	run_mcl.py
    # Input: (1) a path to an MCL input file (2) a path to MCL's output.
    # Output: MCL analysis.
    # Can be parallelized on cluster
    step = '09'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_mcl_analysis'
    script_path = os.path.join(args.src_dir, 'run_mcl.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_run_mcl.txt')
    if not os.path.exists(done_file_path):
        logger.info('Executing MCL...')
        for putative_orthologs_group in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_group_prefix = os.path.splitext(putative_orthologs_group)[0]
            output_file_name = f'{putative_orthologs_group_prefix}.{dir_name}'
            params = [os.path.join(previous_pipeline_step_output_dir, putative_orthologs_group),
                      os.path.join(pipeline_step_output_dir, output_file_name)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules_as_list=['MCL-edge/mcl-14-137'])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 10.	verify_cluster.py
    # Input: (1) mcl analysis file (2) a path to which the file will be moved if relevant (3) optional: maximum number of clusters allowed [default=1]
    # Output: filter irrelevant clusters by moving the relevant to an output directory
    # Can be parallelized on cluster
    step = '10'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_verified_clusters'
    script_path = os.path.join(args.src_dir, 'verify_cluster.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_verify_cluster.txt')
    if not os.path.exists(done_file_path):
        logger.info('Verifying clusters...')
        for putative_orthologs_group in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_group_prefix = os.path.splitext(putative_orthologs_group)[0]
            job_name = os.path.split(putative_orthologs_group_prefix)[-1]
            params = [os.path.join(previous_pipeline_step_output_dir, putative_orthologs_group),
                      os.path.join(pipeline_step_output_dir, putative_orthologs_group)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name, queue_name=args.queue_name)
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 11.	construct_final_orthologs_table.py
    # Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
    # Output: aggregates all the well-clustered OGs to the final table.
    step = '11'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_final_table'
    script_path = os.path.join(args.src_dir, 'construct_final_orthologs_table.py')
    num_of_expected_results = 1 # a single job that prepares all the files
    final_orthologs_table_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    final_orthologs_table_file_path = os.path.join(final_orthologs_table_path, 'final_orthologs_table.csv')
    final_table_header_path = os.path.join(final_orthologs_table_path, 'final_table_header.txt')
    done_file_path = os.path.join(done_files_dir, f'{step}_construct_final_orthologs_table.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing final orthologs table...')
        job_name = os.path.split(final_orthologs_table_file_path)[-1]
        params = [putative_orthologs_table_path,
                  previous_pipeline_step_output_dir,
                  final_orthologs_table_file_path,
                  final_table_header_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name, queue_name=args.queue_name)
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 12.	extract_orthologs_sequences.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    step = '12'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_orthologs_groups_sequences'
    script_path = os.path.join(args.src_dir, 'extract_orthologs_sequences.py')
    num_of_expected_results = 0
    orthologs_sequences_dir_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    num_of_strains_path = os.path.join(final_orthologs_table_path, 'num_of_strains.txt')
    done_file_path = os.path.join(done_files_dir, f'{step}_extract_orthologs_sequences.txt')
    if not os.path.exists(done_file_path):
        logger.info(f'There are {len(os.listdir(previous_pipeline_step_output_dir))} verified clusters in {previous_pipeline_step_output_dir}.')
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')
        logger.debug(f'The verified clusters in {previous_pipeline_step_output_dir} are the following:')
        logger.debug(os.listdir(previous_pipeline_step_output_dir))
        i = 0
        with open(final_orthologs_table_file_path) as f:
            final_table_header = f.readline().rstrip()
            for line in f:
                cluster_members = line.rstrip()
                output_file_name = f'og_{i}.{dir_name}'
                params = [ORFs_dir,
                          f'"{final_table_header}"',  # should be flanked by quotes because it might contain spaces...
                          f'"{cluster_members}"',  # should be flanked by quotes because it might contain spaces...
                          os.path.join(orthologs_sequences_dir_path, output_file_name)]
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
                num_of_expected_results += 1
                i += 1

        # extract number of strains for core genome analysis later on
        num_of_strains = len(final_table_header.split(','))
        with open(num_of_strains_path, 'w') as f:
            f.write(f'{num_of_strains}\n')

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 13.	extract_groups_sizes_frequency
    # Input: TODO
    # Output: TODO
    step = '13'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_groups_sizes_frequency'
    group_sizes_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    groups_sizes_frequency_file_prefix = os.path.join(group_sizes_path, 'groups_sizes_frequency')
    done_file_path = os.path.join(done_files_dir, f'{step}_extract_groups_sizes_frequency.txt')
    if not os.path.exists(done_file_path):
        logger.info('Collecting sizes...')

        group_sizes = []
        with open(final_orthologs_table_file_path) as f:
            final_table_header = f.readline().rstrip()
            for line in f:
                cluster_members = line.rstrip()
                size = sum(bool(item) for item in cluster_members.split(','))  # count non empty entries
                group_sizes.append(str(size))

        groups_sizes_frequency_raw_file_path = groups_sizes_frequency_file_prefix + '.txt'
        with open(groups_sizes_frequency_raw_file_path, 'w') as f:
            f.write('\n'.join(group_sizes)) #f.write('\n'.join([f'{size},{group_size_to_counts_dict[size]}' for size in group_size_to_counts_dict]))

        groups_sizes_frequency_png_file_path = groups_sizes_frequency_file_prefix + '.png'
        generate_bar_plot(groups_sizes_frequency_raw_file_path, groups_sizes_frequency_png_file_path,
            xlabel='\nOrthologs group size', ylabel='Counts\n', dpi=100) #todo: remove dpi

        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 14.	extract_orfs_statistics
    # Input: TODO
    # Output: TODO
    step = '14'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_orfs_statistics'
    orfs_statistics_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    orfs_counts_frequency_file_prefix = os.path.join(orfs_statistics_path, 'orfs_counts')
    orfs_gc_content_file_prefix = os.path.join(orfs_statistics_path, 'orfs_gc_contents')
    done_file_path = os.path.join(done_files_dir, f'{step}_extract_orfs_statistics.txt')
    if not os.path.exists(done_file_path):
        logger.info('Collecting orfs counts...')

        orfs_counts = []
        gc_contents = []
        for file in os.listdir(ORFs_dir):
            orf_path = os.path.join(ORFs_dir, file)
            with open(orf_path) as f:
                num_of_nucleotides = 0
                num_of_GC = 0
                orfs_count = 0
                sequence = ''
                for line in f:
                    if line.startswith('>'):
                        orfs_count += 1
                        num_of_nucleotides += len(sequence)
                        num_of_GC += sequence.count('G') + sequence.count('C')
                    else:
                        sequence += line.rstrip().upper()

                # don't forget last record!!
                orfs_count += 1
                num_of_nucleotides += len(sequence)
                num_of_GC += sequence.count('G') + sequence.count('C')

                orfs_counts.append(str(orfs_count))
                gc_contents.append(str(num_of_GC / num_of_nucleotides))

        orfs_counts_frequency_raw_file_path = orfs_counts_frequency_file_prefix + '.txt'
        with open(orfs_counts_frequency_raw_file_path, 'w') as f:
            f.write('\n'.join(orfs_counts))#[f'{count},{orfs_counts[count]}' for count in orfs_counts]))

        orfs_counts_frequency_png_file_path = orfs_counts_frequency_file_prefix + '.png'
        generate_boxplot(orfs_counts_frequency_raw_file_path, orfs_counts_frequency_png_file_path,
                        xlabel='\nORFs count per genome', dpi=100)

        orfs_gc_content_raw_file_path = orfs_gc_content_file_prefix + '.txt'
        with open(orfs_gc_content_raw_file_path, 'w') as f:
            f.write('\n'.join(gc_contents))

        orfs_gc_content_png_file_path = orfs_gc_content_file_prefix + '.png'
        generate_boxplot(orfs_gc_content_raw_file_path, orfs_gc_content_png_file_path,
                        xlabel='\nGC content per genome', dpi=100)

        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 15.	align_orthologs_group.py
    # Input: TODO
    # Output: TODO
    step = '15'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_aligned_orthologs_groups'
    script_path = os.path.join(args.src_dir, 'align_orthologs_group.py')
    num_of_expected_results = 0
    alignments_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_aligned_orthologs_groups.txt')
    if not os.path.exists(done_file_path):
        logger.info('Aligning orthologs groups...')
        for og_file in os.listdir(orthologs_sequences_dir_path):
            og_path = os.path.join(orthologs_sequences_dir_path, og_file)
            og_file_prefix = os.path.splitext(og_file)[0]
            alignment_path = os.path.join(alignments_path, f'{og_file_prefix}.fasta')
            params = [og_path, alignment_path]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=og_file_prefix,
                                 queue_name=args.queue_name, required_modules_as_list=['mafft/mafft7313'])
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results,
                         error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    #16.	extract aligned core genome.py
    step = '16'
    logger.info(f'Step {step}: {"_" * 100}')
    dir_name = f'{step}_aligned_core_genome'
    script_path = os.path.join(args.src_dir, 'extract_core_genome.py')
    num_of_expected_results = 1
    aligned_core_genome_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_extract_aligned_core_genome.txt')
    aligned_core_genome_file_path = os.path.join(aligned_core_genome_path, 'aligned_core_genome.fasta')
    if not os.path.exists(done_file_path):
        logger.info('Extracting aligned core genome...')
        with open(num_of_strains_path) as f:
            num_of_strains = f.read().rstrip()

        params = [alignments_path,
                  num_of_strains,
                  aligned_core_genome_file_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name='core_genome',
                             queue_name=args.queue_name)

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results,
                         error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 17.	reconstruct_species_phylogeny.py
    step = '17'
    logger.info(f'Step {step}: {"_" * 100}')
    dir_name = f'{step}_species_phylogeny'
    script_path = os.path.join(args.src_dir, 'reconstruct_species_phylogeny.py')
    num_of_expected_results = 1
    phylogeny_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_species_phylogeny.txt')
    if not os.path.exists(done_file_path):
        logger.info('Reconstructing species phylogeny...')

        output_file_prefix = 'species_tree'
        phylogenetic_raw_tree_path = os.path.join(phylogeny_path, output_file_prefix + '.txt')
        params = [aligned_core_genome_file_path,
                  phylogenetic_raw_tree_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_prefix,
                             queue_name=args.queue_name, required_modules_as_list=['RAxML/RAxML-8.1.3'])

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results,
                         error_file_path)

        phylogenetic_png_tree_path = os.path.join(phylogeny_path, output_file_prefix + '.png')
        generate_tree_plot(phylogenetic_raw_tree_path, phylogenetic_png_tree_path)

        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # Final step: gather relevant results, zip them together and update html file
    logger.info(f'FINAL STEP: {"_"*100}')
    final_output_dir_name = f'{CONSTS.WEBSERVER_NAME}_outputs'
    final_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, final_output_dir_name)
    done_file_path = os.path.join(done_files_dir, final_output_dir_name + '.txt')
    if not os.path.exists(done_file_path):
        logger.info('Gathering results to final output dir...')

        # move orthologs table
        add_results_to_final_dir(final_orthologs_table_path, final_output_dir)

        # move unaligned sequences
        add_results_to_final_dir(orthologs_sequences_dir_path, final_output_dir)

        # move groups sizes
        add_results_to_final_dir(group_sizes_path, final_output_dir)

        # move orfs statistics
        add_results_to_final_dir(orfs_statistics_path, final_output_dir)

        # move aligned sequences
        add_results_to_final_dir(alignments_path, final_output_dir)

        # move core genome dir
        add_results_to_final_dir(aligned_core_genome_path, final_output_dir)

        # move species tree dir
        add_results_to_final_dir(phylogeny_path, final_output_dir)

        logger.info('Zipping results folder...')
        shutil.make_archive(final_output_dir, 'zip', final_output_dir)

        logger.info(f'Moving results to parent dir... ({meta_output_dir})')
        try:
            shutil.move(final_output_dir+'.zip', meta_output_dir)
        except shutil.Error as e:
            logger.error(e.args[0])
        try:
            shutil.move(final_output_dir, meta_output_dir)
        except shutil.Error as e:
            logger.error(e.args[0])

        #TODO: remove HERE rest of intermediate outputs

        file_writer.write_to_file(done_file_path)
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')

    logger.info('Editing results html...')
    edit_success_html(output_html_path, meta_output_dir, final_output_dir_name, run_number, CONSTS)

    status = 'is done'


except Exception as e:
    status = 'was failed'
    from time import ctime
    import os
    import logging
    logger = logging.getLogger('main')  # use logger instead of printing

    msg = 'M1CR0B1AL1Z3R failed :('

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    logger.error(f'\n\n{"$" * 100}\n\n{msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\ne.args[0]: {e.args[0]}\n\n{"$" * 100}')

    edit_failure_html(output_html_path, run_number, msg, CONSTS)

    notify_admin(meta_output_dir, meta_output_url, run_number, CONSTS)

end = time()

results_location = output_url if remote_run else args.output_dir
msg = f'M1CR0B1AL1Z3R pipeline {status}'
if status == 'is done':
    msg += f' (Took {measure_time(int(end-start))}).\nResults can be found at {results_location}.'
    if remote_run and run_number.lower() != 'example' and False: #TODO: remove the "and False" once ready.
        # remove raw data from the server
        try:
            shutil.rmtree(data_path)
        except:
            pass
        # remove raw data from the server
        try:
            shutil.rmtree(args.output_dir)
        except:
            pass
else:
    msg += f'.\nFor further information please visit: {results_location}'
logger.info(msg)
send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>', args.email, subject=f'Microbialzer {status}.', content=msg)

