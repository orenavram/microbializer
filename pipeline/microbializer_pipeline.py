try:
    import argparse
    import sys
    import tarfile
    from sh import gunzip

    sys.path.insert(0, "/groups/pupko/orenavr2/microbializer/") #ADD auxiliaries path
    sys.path.insert(0, "/groups/pupko/orenavr2/microbializer/auxiliaries/") #ADD file_writer
    from auxiliaries import file_writer #ADD file_writer
    from auxiliaries.email_sender import send_email
    from auxiliaries.pipeline_auxiliaries import *

    start = time()

    parser = argparse.ArgumentParser()
    parser.add_argument('contigs_dir', help='path to a folder with the genomic sequences. This folder may be zipped, as well the files in it.',
                        type=lambda file_path: str(file_path) if os.path.exists(file_path) else parser.error(f'{file_path} does not exist!'))
    parser.add_argument('output_dir', help='directory where the output files will be written to')
    parser.add_argument('email', help='A notification will be sent once the pipeline is done',
                        default='orenavram@gmail.com')
    parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                        choices=['pupko', 'itaym', 'lilach', 'bioseq', 'bental'], default='pupko')
    parser.add_argument('--dummy_delimiter',
                        help='The queue does not "like" very long commands. A dummy delimiter is used to break each row into different commands of a single job',
                        default='!@#')
    parser.add_argument('--src_dir', help='source code directory', type=lambda s: s.rstrip('/'), default='/groups/pupko/orenavr2/microbializer/pipeline/')
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
    error_file_path = os.path.join(args.output_dir, 'error.txt')

    tmp_dir = os.path.join(args.output_dir, 'tmp_dir')
    create_dir(tmp_dir)

    done_files_dir = os.path.join(args.output_dir, 'done')
    create_dir(done_files_dir)

    data_path = args .contigs_dir
    logger.info(data_path)

    # extract tar folder
    if not os.path.isdir(data_path) and tarfile.is_tarfile(data_path):
        with tarfile.open(data_path, 'r:gz') as f:
            f.extractall() #unzip tar folder
        data_path = data_path.split('.tar')[0] # e.g., /groups/pupko/orenavr2/microbializer/example_data.tar.gz
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
    logger.warning(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_ORFs'
    script_path = os.path.join(args.src_dir, 'extract_orfs.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_extract_orfs.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Extracting ORFs...')
        for fasta_file in os.listdir(data_path):
            fasta_file_prefix = os.path.splitext(fasta_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            output_coord_name = f'{fasta_file_prefix}.gene_coordinates'
            params = [os.path.join(data_path, fasta_file), os.path.join(pipeline_step_output_dir, output_file_name),
                      os.path.join(pipeline_step_output_dir, output_coord_name)] #Shir - path to translated sequences file
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules_as_list=['prodigal/prodigal-2.6.3'])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 2.  create_blast_DB.py
    # Input: path to gene file to create DB from
    # Output: Blast DB of the input file
    # Can be parallelized on cluster
    step = '02'
    logger.warning(f'Step {step}: {"_"*100}')
    ORFs_dir = previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_blast_db'
    script_path = os.path.join(args.src_dir, 'create_blast_DB.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_create_blast_DB.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Creating Blast DB...')
        for fasta_file in os.listdir(previous_pipeline_step_output_dir):
            file_path = os.path.join(previous_pipeline_step_output_dir, fasta_file)
            fasta_file_prefix = os.path.splitext(fasta_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            params = [file_path, os.path.join(pipeline_step_output_dir,output_file_name)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=fasta_file_prefix,queue_name=args.queue_name, required_modules_as_list=['blast/blast-2.2.30'])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 3.	blast_all_vs_all.py
    # Input: (1) 2 input paths for 2 (different) genes files, g1 and g2 (2) an output file path (with a suffix as follows: i_vs_j_blast.tsv. especially relevant for the wrapper).
    # Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
    # Precisely, for each gene x in g1, blast x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
    # Can be parallelized on cluster
    step = '03'
    logger.warning(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_blast_analysis'
    script_path = os.path.join(args.src_dir, 'blast_all_vs_all.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_blast_all_vs_all.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Blasting...')
        blasted_pairs = set()
        for fasta_file_name in os.listdir(ORFs_dir):
            fasta_file_path = os.path.join(ORFs_dir, fasta_file_name)
            strain1_name = os.path.splitext(fasta_file_name)[0]
            for db_file in os.listdir(previous_pipeline_step_output_dir):
                # 2 extensions should be removed GCF_000006945.2_ASM694v2_cds_from_genomic.02_blast_db.nhr or GCF_000006945.2_ASM694v2_cds_from_genomic.02_blast_db.nin
                strain2_name = os.path.splitext(os.path.splitext(db_file)[0])[0]
                pair = (strain1_name, strain2_name)
                if pair not in blasted_pairs and strain1_name != strain2_name:
                    blasted_pairs.add(pair)
                    logger.debug(f'{"#"*100}\nBlasting pair: {strain1_name}, {strain2_name}')
                    logger.debug(blasted_pairs)
                    db_prefix = os.path.splitext(os.path.join(previous_pipeline_step_output_dir, db_file))[0]
                    output_file_name = f'{strain1_name}_vs_{strain2_name}.{dir_name}'
                    params = [fasta_file_path, db_prefix, os.path.join(pipeline_step_output_dir, output_file_name)]
                    submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules_as_list=['blast/blast-2.2.30'])
                    num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 4.	filter_blast.py
    # Input: (1) a path for a i_vs_j_blast.tsv file (2) an output path (with a suffix as follows: i_vs_j_filtered.tsv. especially relevant for the wrapper).
    # Output: the same format of the input file containing only pairs that passed the filtration. For each row in the input file (pair of genes), apply the following filters:
    # 1. at least X% similarity
    # 2. at least X% of the length
    # 3.# write each pair to the output file if it passed all the above filters.
    # Can be parallelized on cluster
    step = '04'
    logger.warning(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_blast_filtered'
    script_path = os.path.join(args.src_dir, 'filter_blast_results.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_filter_blast_results.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Filtering...')
        for blast_results_file in os.listdir(previous_pipeline_step_output_dir):
            fasta_file_prefix = os.path.splitext(blast_results_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            params = [os.path.join(previous_pipeline_step_output_dir, blast_results_file), os.path.join(pipeline_step_output_dir, output_file_name)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 5.	find_reciprocal_hits.py
    # Input: (1) a path for a i_vs_j.blast_filtered file and a path for a j_vs_i.blast_filtered file (2) an output path (with a suffix as follows: i_vs_j_reciprocal_hits.tsv
    # Output: a tab delimited file containing only reciprocal best-hit pairs and their bit score from blast, i.e., if x’s best hit was y in the first file with bit score z, x\ty\tz will appear in the output file only if y’s best hit in the other file will be x.
    # Can be parallelized on cluster
    step = '05'
    logger.warning(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_reciprocal_hits'
    script_path = os.path.join(args.src_dir, 'find_reciprocal_hits.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_find_reciprocal_hits.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Extracting reciprocal hits...')
        for blast_filtered_results_file in os.listdir(previous_pipeline_step_output_dir):
            file_name_prefix = os.path.splitext(blast_filtered_results_file)[0]
            strain1, strain2 = file_name_prefix.split('_vs_')[:2]
            logger.debug(f'\nstrain1: {strain1}\nstrain2: {strain2}')
            if strain1<strain2:
                # avoid both id1,id2 and id2,id1
                reciprocal_file_name_prefix = "_vs_".join([strain2, strain1])
                reciprocal_file_name = reciprocal_file_name_prefix + os.path.splitext(blast_filtered_results_file)[1]
                output_file_name = f'{file_name_prefix}.{dir_name}'
                params = [os.path.join(previous_pipeline_step_output_dir, blast_filtered_results_file),
                          os.path.join(previous_pipeline_step_output_dir, reciprocal_file_name),
                          os.path.join(pipeline_step_output_dir, output_file_name)]
                #print("\n\n\n",script_path,params[0], params[1] , params[2] ,"\n\n\n\n")
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
                num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 6. concatenate_reciprocal_hits
    # Input: path to folder with all reciprocal hits files
    # Output: concatenated file of all reciprocal hits files
    # CANNOT be parallelized on cluster
    step = '06'
    logger.warning(f'Step {step}: {"_"*100}')
    all_reciprocal_hits_file = os.path.join(args.output_dir, 'concatenated_all_reciprocal_hits.txt')
    done_file_path = os.path.join(done_files_dir, f'{step}_concatenate_reciprocal_hits.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Concatenating reciprocal hits...')
        cmd = f'cat {pipeline_step_output_dir}/*.{dir_name} > {all_reciprocal_hits_file}'
        subprocess.run(cmd, shell = True)
        # No need to wait...
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 7.	construct_putative_orthologs_table.py
    # Input: (1) a path for a i_vs_j_reciprocal_hits.tsv file (2) a path for a putative orthologs file (with a single line).
    # Output: updates the table with the info from the reciprocal hit file.
    # CANNOT be parallelized on cluster
    step = '07'
    logger.warning(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_putative_table'
    script_path = os.path.join(args.src_dir, 'construct_putative_orthologs_table.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    putative_orthologs_table_path = os.path.join(pipeline_step_output_dir, 'putative_orthologs_table.txt')
    done_file_path = os.path.join(done_files_dir, f'{step}_construct_putative_orthologs_table.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Constructing putative orthologs table...')
        job_name = os.path.split(script_path)[-1]
        params = [all_reciprocal_hits_file, putative_orthologs_table_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name,
                             queue_name=args.queue_name)
        num_of_expected_results = 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # # 6.	split_putative_orthologs_table
    # # Input: (1) a path for a putative orthologs table (each row in this table is a putative orthologous set) (2) an output_path to a directory.
    # # Output: each line is written to a separate file in the output directory.
    # # Can be parallelized on cluster (better as a subprocess)
    # logger.warning(f'Step 6: {"_"*100}')
    # dir_name = 'splitted_putative_orthologs_table'
    # # script_path = os.path.join(args.src_dir, '..', 'auxiliaries', 'file_writer.py')
    # # num_of_expected_results = 0
    # pipeline_step_output_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)[0]  # pipeline_step_tmp_dir is not needed
    # done_file_path = os.path.join(done_files_dir, '6_split_putative_orthologs_table.txt')
    # if not os.path.exists(done_file_path):
    #     logger.warning('Splitting putative orthologs table...')
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
    #     logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 8   prepare_files_for_mcl.py
    # Input: (1) a path for a concatenated all reciprocal hits file (2) a path for a putative orthologs file (3) a path for an output folder
    # Output: an input file for MCL for each putative orthologs group
    # CANNOT be parallelized on cluster (if running on the concatenated file)
    step = '08'
    logger.warning(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_mcl_input_files'
    script_path = os.path.join(args.src_dir, 'prepare_files_for_mcl.py')
    num_of_expected_results = 1 # a single job that prepares all the files
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_prepare_files_for_mcl.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Preparing files for MCL...')
        job_name = os.path.split(script_path)[-1]
        params = [all_reciprocal_hits_file, putative_orthologs_table_path, pipeline_step_output_dir]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name, queue_name=args.queue_name)
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 9.	run_mcl.py
    # Input: (1) a path to an MCL input file (2) a path to MCL's output.
    # Output: MCL analysis.
    # Can be parallelized on cluster
    step = '09'
    logger.warning(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_mcl_analysis'
    script_path = os.path.join(args.src_dir, 'run_mcl.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_run_mcl.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Executing MCL...')
        for putative_orthologs_set in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_set_prefix = os.path.splitext(putative_orthologs_set)[0]
            output_file_name = f'{putative_orthologs_set_prefix}.{dir_name}'
            params = [os.path.join(previous_pipeline_step_output_dir, putative_orthologs_set),
                      os.path.join(pipeline_step_output_dir, output_file_name)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name, required_modules_as_list=['MCL-edge/mcl-14-137'])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 10.	verify_cluster.py
    # Input: (1) mcl analysis file (2) a path to which the file will be moved if relevant (3) optional: maximum number of clusters allowed [default=1]
    # Output: filter irrelevant clusters by moving the relevant to an output directory
    # Can be parallelized on cluster
    step = '10'
    logger.warning(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_verified_clusters'
    script_path = os.path.join(args.src_dir, 'verify_cluster.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_verify_cluster.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Verifying clusters...')
        for putative_orthologs_set in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_set_prefix = os.path.splitext(putative_orthologs_set)[0]
            job_name = os.path.split(putative_orthologs_set_prefix)[-1]
            params = [os.path.join(previous_pipeline_step_output_dir, putative_orthologs_set),
                      os.path.join(pipeline_step_output_dir, putative_orthologs_set)]
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name, queue_name=args.queue_name, required_modules_as_list=['MCL-edge/mcl-14-137'])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 11.	construct_final_orthologs_table.py
    # Input: (1) a path for directory with all the verified OGs (2) an output path to a final OGs table.
    # Output: aggregates all the well-clustered OGs to the final table.
    step = '11'
    logger.warning(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_final_table'
    script_path = os.path.join(args.src_dir, 'construct_final_orthologs_table.py')
    num_of_expected_results = 1 # a single job that prepares all the files
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    final_orthologs_table_path = os.path.join(pipeline_step_output_dir, 'final_orthologs_table.txt')
    final_table_header_path = os.path.join(pipeline_step_output_dir, 'final_table_header.txt')
    done_file_path = os.path.join(done_files_dir, f'{step}_construct_final_orthologs_table.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Constructing final orthologs table...')
        job_name = os.path.split(final_orthologs_table_path)[-1]
        params = [putative_orthologs_table_path,
                  previous_pipeline_step_output_dir,
                  final_orthologs_table_path,
                  final_table_header_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name, queue_name=args.queue_name)
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')


    # 12.	extract_orthologs_sequences.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    step = '12'
    logger.warning(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_orthologs_sets_sequences'
    script_path = os.path.join(args.src_dir, 'extract_orthologs_sequences.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{step}_extract_orthologs_sequences.txt')
    if not os.path.exists(done_file_path):
        logger.warning('Extracting orthologs set sequences according to final orthologs table...')
        logger.warning(f'The verified clusters in {previous_pipeline_step_output_dir} are the following:')
        logger.warning(os.listdir(previous_pipeline_step_output_dir))
        # with open(final_table_header_path) as f:
        #     final_table_header = f.readline().rstrip()

        i = 0
        with open(final_orthologs_table_path) as f:
            final_table_header = f.readline().rstrip()
            for line in f:
                cluster_members = line.rstrip()
                output_file_name = f'og_{i}.{dir_name}'
                params = [ORFs_dir,
                          f'"{final_table_header}"',  # should be flanked by quotes because it might contain spaces...
                          f'"{cluster_members}"',  # should be flanked by quotes because it might contain spaces...
                          os.path.join(pipeline_step_output_dir, output_file_name)]
                print(f'params for orthologs_sets_sequences: {params}')
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name, queue_name=args.queue_name)
                num_of_expected_results += 1
                i += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path)
    else:
        logger.warning(f'done file {done_file_path} already exists. Skipping step...')

    status = 'is done'

except Exception as e:
    from time import ctime
    import logging
    logger = logging.getLogger('main')  # use logger instead of printing

    status = 'was failed'
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    logger.error('\n\n' + '$' * 60 + '\n\n')
    logger.error(f'{ctime()}: microbializer failed :(\n\n')
    logger.error(f'{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\n')
    logger.error('$' * 60)

    raise e
    # TODO: print message from error_file_path
    # with open(error_file_path, 'w') as f:
    #     f.write(e.args[0])

end = time()
msg = f'microbializer pipeline {status}. Results can be found at {args.output_dir}. Took {measure_time(int(end-start))}'
logger.warning(msg)
send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>', args.email, subject=f'Microbialzer {status}.', content=msg)

