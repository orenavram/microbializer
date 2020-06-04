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
    file_with_job_id_on_qstat = os.path.join(meta_output_dir, 'qsub.log')
    if os.path.exists(file_with_job_id_on_qstat):
        with open(file_with_job_id_on_qstat) as f:
            job_id_on_qstat = f.read().strip()
        error_log_path = os.path.join(meta_output_dir, f'{job_id_on_qstat}.ER')
        # TODO: change to ER url and add reading permissions
    # Send me a notification email every time there's a failure
    send_email(smtp_server=CONSTS.SMTP_SERVER,
               sender=CONSTS.ADMIN_EMAIL,
               receiver=CONSTS.OWNER_EMAIL,
               subject=f'{CONSTS.WEBSERVER_NAME} job {run_number} by {email} has been failed: ',
               content=f"{email}\n\n{os.path.join(meta_output_url, 'output.html')}\n\n"
               f"{os.path.join(meta_output_url, 'cgi_debug.txt')}\n\n"
               f"{os.path.join(meta_output_url, error_log_path)}\n\n"
               f"{os.path.join(meta_output_dir, error_log_path.replace('ER', 'OU'))}")


def verify_fasta_format(data_path):
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() + Bio.SeqUtils.IUPACData.ambiguous_dna_letters)
    for file_name in os.listdir(data_path):
        if file_name.endswith('zip') or file_name.endswith('gz') or file_name.endswith('rar'):
            return f'{file_name} is a binary file (rather than textual). Please upload your genomes as text files (such as fas or fna).'
        file_path = os.path.join(data_path, file_name)
        strain_name = os.path.splitext(file_name)[0]
        with open(file_path) as f:
            line_number = 0
            try:
                line = f.readline()
                line_number += 1
                if not line.startswith('>'):
                    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. First line in "{file_name}" starts with "{line[0]}" instead of ">".'
                previous_line_was_header = True
                putative_end_of_file = False
                curated_content = f'>{strain_name}_{line[1:]}'.replace("|", "_")
                for line in f:
                    line_number += 1
                    line = line.strip()
                    if not line:
                        if not putative_end_of_file: # ignore trailing empty lines
                            putative_end_of_file = line_number
                        continue
                    if putative_end_of_file:  # non empty line after empty line
                        return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {putative_end_of_file} in "{file_name}" is empty.'
                    if line.startswith('>'):
                        if previous_line_was_header:
                            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. "{file_name}" contains an empty record. Both lines {line_number-1} and {line_number} start with ">".'
                        else:
                            previous_line_was_header = True
                            curated_content += f'>{strain_name}_{line[1:]}\n'.replace("|", "_")
                            continue
                    else:  # not a header
                        previous_line_was_header = False
                        for c in line:
                            if c not in legal_chars:
                                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in "{file_name}" contains illegal DNA character "{c}".'
                        curated_content += f'{line}\n'
            except UnicodeDecodeError as e:
                logger.info(e.args)
                line_number += 1  # the line that was failed to be read
                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in "{file_name}" contains one (or more) non <a href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s).'
        # override the old file with the curated content
        with open(file_path, 'w') as f:
            f.write(curated_content)


def edit_progress(output_html_path, progress=None, active=True):
    result = ''
    with open(output_html_path) as f:
        for line in f:
            if 'progress-bar' in line:
                if progress:
                    line = line.split('style')[0]  # <div class="progress-bar ... style="width:0%">\n
                    line += f'style="width:{progress}%">\n'
                if not active:
                    line = line.replace(' active', '')  # <div class="progress-bar progress-bar-striped bg-success active" ...
            result += line

    with open(output_html_path, 'w') as f:
        f.write(result)

def add_results_to_final_dir(source, final_output_dir, copy=True):
    dest = os.path.join(final_output_dir, os.path.split(source)[1])

    try:
        if not copy:
            logger.info(f'Moving {source} TO {dest}')
            shutil.move(source, dest)
        else:
            logger.info(f'Copying {source} TO {dest}')
            shutil.copytree(source, dest)
    except FileExistsError:
        pass

    return dest


def remove_path(path_to_remove):
    logger.info(f'Removing {path_to_remove} ...')
    try:
        shutil.rmtree(path_to_remove)  # maybe it's a folder
    except:
        pass
    try:
        os.remove(path_to_remove)
    except:
        pass


def remove_bootstrap_values(in_tree_path, out_tree_path):
    with open(in_tree_path) as f:
        tree_as_str = f.read()
    import re
    tree_as_str = re.sub('\)\d+:', '):', tree_as_str)
    with open(out_tree_path, 'w') as f:
        f.write(tree_as_str)


def unpack_data(data_path):
    if not os.path.isdir(data_path):
        unzipped_data_path = os.path.join(meta_output_dir, 'data')
        try:
            if tarfile.is_tarfile(data_path):
                logger.info('UnTARing')
                with tarfile.open(data_path, 'r:gz') as f:
                    f.extractall(path=unzipped_data_path)  # unzip tar folder to parent dir
                logger.info('Succeeded!')
                # data_path = data_path.split('.tar')[0] # e.g., /groups/pupko/orenavr2/microbializer/example_data.tar.gz
                # logger.info(f'Updated data_path is:\n{data_path}')
            elif data_path.endswith('.gz'):  # gunzip gz file
                subprocess.run(f'gunzip -f "{data_path}"', shell=True)
                unzipped_data_path = data_path[:-3]  # trim the ".gz"
            else:
                logger.info('UnZIPing')
                shutil.unpack_archive(data_path, extract_dir=unzipped_data_path)  # unzip tar folder to parent dir
        except Exception as e:
            logger.info(e)
            remove_path(data_path)
            fail(f'Illegal file format. Please upload either a '
                 f'<a href="https://support.microsoft.com/en-us/help/14200/windows-compress-uncompress-zip-files" target="_blank">.zip</a> file or a '
                 f'<a href="https://linhost.info/2012/08/gzip-files-in-windows/" target="_blank">.tar.gz</a> file in which each file is a '
                 f'<a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a> containing genomic sequence of a different species',
                 error_file_path)
        logger.info('Succeeded!')
        # data_path = os.path.join(meta_output_dir, 'data') # e.g., /groups/pupko/orenavr2/microbializer/example_data.tar.gz
        # logger.info(f'Updated data_path is:\n{data_path}')

        if not os.path.exists(unzipped_data_path):
            fail(f'Failed to unzip {os.path.split(data_path)[-1]} (maybe it is empty?)', error_file_path)

        if not os.path.isdir(unzipped_data_path):
            fail('Archived file content is not a folder', error_file_path)

        file = [x for x in os.listdir(unzipped_data_path) if not x.startswith(('_', '.'))][0]
        logger.info(f'first file in {unzipped_data_path} is:\n{file}')
        if os.path.isdir(os.path.join(unzipped_data_path, file)):
            data_path = os.path.join(unzipped_data_path, file)
            file = [x for x in os.listdir(data_path) if not x.startswith(('_', '.'))][0]
            if os.path.isdir(os.path.join(data_path, file)):
                fail('More than a 2-levels folder...', error_file_path)
        else:
            data_path = unzipped_data_path

    logger.info(f'Updated data_path is:\n{data_path}')
    for file in os.listdir(data_path):
        file_path = os.path.join(data_path, file)
        if file_path.endswith('gz'):  # gunzip gz files in $data_path if any
            subprocess.run(f'gunzip -f "{file_path}"', shell=True)

    return data_path


try:
    import argparse
    import sys
    import os
    import tarfile
    import shutil
    import mmap
    import Bio.SeqUtils

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

    import CONSTANTS as CONSTS

    from html_editor import edit_success_html, edit_failure_html

    start = time()

    parser = argparse.ArgumentParser()
    parser.add_argument('contigs_dir', help='path to a folder with the genomic sequences. This folder may be zipped, as well the files in it.',
                        type=lambda path: path.rstrip('/') if os.path.exists(path) else parser.error(f'{path} does not exist!'))
    parser.add_argument('output_dir', help='directory where the output files will be written to',
                        type=lambda path: path.rstrip('/'))
    parser.add_argument('--email', help='A notification will be sent once the pipeline is done',
                        default=CONSTS.OWNER_EMAIL)
    parser.add_argument('--identity_cutoff', default=80, type=lambda x: eval(x),
                        help='minimum required percent of identity level (lower values will be filtered out)')
    parser.add_argument('--e_value_cutoff', default=0.01, type=lambda x: eval(x), # if 0 <= eval(x) <= 100 else parser.error(f"Can't use {x} as percent!"),
                        help='maxmimum permitted e-value (0 <= e_value_cutoff <= 1; higher values will be filtered out).')
    parser.add_argument('--core_minimal_percentage', default=99.9, type=lambda x: eval(x), # if 0 <= eval(x) <= 100 else parser.error(f"Can't use {x} as percent!"),
                        help='the minimum required percent of gene members that is needed to be considered a core gene. For example: (1) 100 means that for a gene to be considered core, all strains should have a member in the group.\n(2) 50 means that for a gene to be considered core, at least half of the strains should have a member in the group.\n(3) 0 means that every gene should be considered as a core gene.')
    parser.add_argument('--bootstrap', default='no', choices=['yes', 'no'],
                        help='whether or not to apply bootstrap procedure over the reconstructed species tree.')
    parser.add_argument('--outgroup', default=None,
                        help='whether or not to root the species phylogeny.')
    parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to', default='pupkotmpr')  # , choices=['pupkoweb', 'pupkowebr', 'pupkolab', 'pupkolabr', 'pupkotmp', 'pupkotmpr', 'itaym', 'lilach', 'bioseq', 'bental', 'oren.q', 'bioseq20.q'])
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

    meta_output_dir = os.path.join(os.path.split(args.output_dir)[0])
    logger.info(f'meta_output_dir is: {meta_output_dir}')

    error_file_path = os.path.join(meta_output_dir, 'error.txt')

    run_number = os.path.join(os.path.split(meta_output_dir)[1])
    logger.info(f'run_number is {run_number}')

    output_html_path = os.path.join(meta_output_dir, 'output.html')
    logger.info(f'output_html_path is {output_html_path}')

    output_url = os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number, 'output.html')
    logger.info(f'output_url is {output_url}')

    meta_output_url = os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number)

    tmp_dir = os.path.join(args.output_dir, 'tmp_dir')
    create_dir(tmp_dir)

    done_files_dir = os.path.join(args.output_dir, 'done')
    create_dir(done_files_dir)

    data_path = args.contigs_dir
    logger.info(f'data_path is: {data_path}')

    delimiter = ','

    # extract zip and detect data folder
    data_path = unpack_data(data_path)

    for system_file in os.listdir(data_path):
        if system_file.startswith(('.', '_')):
            system_file_path = os.path.join(data_path, system_file)
            logger.warning(f'Removing system file: {system_file_path}')
            remove_path(system_file_path)

    # have to be AFTER system files removal (in the weird case a file name starts with a space)
    filename_prefixes = set()
    illegal_chars = ' ();'
    for file_name in os.listdir(data_path):

        filename_prefix = os.path.splitext(file_name)[0]
        if filename_prefix in filename_prefixes:
            error_msg = f'Two (or more) of the uploaded geonmes contain the same name (prefix), e.g., {filename_prefix}. Please make sure each file name is unique.'
            fail(error_msg, error_file_path)
        filename_prefixes.add(filename_prefix)

        for char in illegal_chars:
            if char in file_name:
                new_file_name = file_name.replace(char, '_')
                logger.info(f'File name with illegal character "{char}" was detected!\n'
                            f'Renaming file path from:\n'
                            f'{os.path.join(data_path, file_name)}\n'
                            f'to this:\n'
                            f'{os.path.join(data_path, new_file_name)}')
                try:
                    os.rename(os.path.join(data_path, file_name), os.path.join(data_path, new_file_name))
                    file_name = new_file_name
                except:
                    error_msg = f'One (or more) file name(s) contain illegal character "{char}".<br>\nIn order to avoid downstream parsing errors, {CONSTS.WEBSERVER_NAME} automatically replaces these spaces with dashes. For some reason, the replacement for {file_name} failed. Please remove space characters from your file names and re-submit your job.'
                    fail(error_msg, error_file_path)

    number_of_genomes = len(os.listdir(data_path))
    logger.info(f'Number of genomes to analyze is {number_of_genomes}')
    logger.info(f'data_path contains the following {number_of_genomes} files:\n{os.listdir(data_path)}')

    # check MINimal number of genomes
    min_number_of_genomes_to_analyze = 2
    if number_of_genomes < min_number_of_genomes_to_analyze:
        error_msg = f'The dataset contains too few genomes ({CONSTS.WEBSERVER_NAME} does comparative analysis and thus needs at least 2 genomes).'
        fail(error_msg, error_file_path)

    # check MAXimal number of genomes
    max_number_of_genomes_to_analyze = 250
    if number_of_genomes > max_number_of_genomes_to_analyze and args.outgroup.lower() != 'michael':
        error_msg = f'The dataset contains too many genomes (currently {CONSTS.WEBSERVER_NAME} is able to handle up to 200 genomes due to technical issues. We are working to optimize it and soon it will be able to handle more genomes!).'
        fail(error_msg, error_file_path)

    # must be only after the spaces removal from the species names!!
    verification_error = verify_fasta_format(data_path)
    if verification_error:
        remove_path(data_path)
        fail(verification_error, error_file_path)

    # 1.	extract_orfs_sequences.py
    # Input: (1) an input path for a fasta file with contigs/full genome (2) an output file path (with a suffix as follows: i_genes.fasta. especially relevant for the wrapper).
    # Output: a fasta file where under each header, there’s a single gene.
    # Can be parallelized on cluster
    # Prodigal ignores newlines in the middle of a sequence, namely, >bac1\nAAA\nAA\nTTT >bac1\nAAAAATTT will be analyzed identically.
    step = '01'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_ORFs'
    script_path = os.path.join(args.src_dir, 'extract_orfs_sequences.py')
    num_of_expected_results = 0
    ORFs_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Extracting ORFs...')
        for fasta_file in os.listdir(data_path):
            fasta_file_prefix = os.path.splitext(fasta_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            output_coord_name = f'{fasta_file_prefix}.gene_coordinates'
            params = [f'"{os.path.join(data_path, fasta_file)}"',
                      os.path.join(ORFs_dir, output_file_name),
                      os.path.join(ORFs_dir, output_coord_name)] #Shir - path to translated sequences file
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name,
                                 queue_name=args.queue_name, required_modules_as_list=[CONSTS.PRODIGAL])
            num_of_expected_results += 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=10)


    # make sure that all ORF files contain something....
    for file in os.listdir(ORFs_dir):
        try:
            with open(os.path.join(ORFs_dir, file), 'rb', 0) as orf_f, mmap.mmap(orf_f.fileno(), 0, access=mmap.ACCESS_READ) as s:
                if s.find(b'>') > -1:
                    continue
        except:
            error_msg = f'{CONSTS.WEBSERVER_NAME} could not find any ORFs in some of the genomes you provided (e.g., {os.path.splitext(file)[0]}). ' \
                        f'It is recommended to use files that contain at least 20K base pairs (each).'
            fail(error_msg, error_file_path)


    # 2.  create_mmseqs2_DB.py
    # Input: path to gene file to create DB from
    # Output: DB of the input file for mmseqs2
    # Can be parallelized on cluster
    step = '02'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = ORFs_dir
    dir_name = f'{step}_dbs'
    script_path = os.path.join(args.src_dir, 'create_mmseqs2_DB.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 2
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info('Creating DBs...')
        for fasta_file in os.listdir(previous_pipeline_step_output_dir):
            file_path = os.path.join(previous_pipeline_step_output_dir, fasta_file)
            fasta_file_prefix = os.path.splitext(fasta_file)[0]
            output_file_name = f'{fasta_file_prefix}'
            if num_of_aggregated_params > 0:  # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)

            output_prefix = os.path.join(pipeline_step_output_dir, output_file_name)
            params = [file_path,
                      output_prefix,
                      output_prefix, #instead of tmp_dir
                      '-t'] #translate to peptides # Should we let the user decide?

            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=fasta_file_prefix,
                                     queue_name=args.queue_name, more_cmds=more_cmds, required_modules_as_list=[CONSTS.MMSEQS])
                num_of_expected_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=fasta_file_prefix,
                                 queue_name=args.queue_name, more_cmds=more_cmds, required_modules_as_list=[CONSTS.MMSEQS])
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=15)

    # send a subjob that removes all mmseqs intermediate *FOLDERS* (e.g., 3465136234521948 etc..) in tmp_dir
    # does not remove (sge/cmds/log) files. only folders.
    submit_pipeline_step(os.path.join(args.src_dir, 'remove_tmp_folders.py'), [pipeline_step_tmp_dir], pipeline_step_tmp_dir,
                         job_name='remove_dirs_from_tmp', queue_name=args.queue_name)


    # 3.	mmseqs2_all_vs_all.py
    # Input: (1) 2 input paths for 2 (different) genome files (query and target), g1 and g2
    #        (2) an output file path (with a suffix as follows: i_vs_j.tsv. especially relevant for the wrapper).
    # Output: a tsv file where the first column contains g1 genes and the second column includes the corresponding best match.
    # Precisely, for each gene x in g1, query x among all the genes of g2. Let y be the gene in g2 that is the most similar to x among all g2 genes. Append a row to the output file with: ‘{x}\t{y}’.
    # Can be parallelized on cluster
    step = '03'
    logger.info(f'Step {step}: {"_"*100}')
    previous_pipeline_step_output_dir = pipeline_step_output_dir
    dir_name = f'{step}_all_vs_all_analysis'
    script_path = os.path.join(args.src_dir, 'mmseqs2_all_vs_all.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 100 if len(os.listdir(ORFs_dir)) > 100 else 50
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info(f'Querying all VS all (using mmseqs2)...')
        for db1 in os.listdir(ORFs_dir):
            fasta_file_path = os.path.join(ORFs_dir, db1)
            strain1_name = os.path.splitext(db1)[0]
            for db2 in os.listdir(ORFs_dir):
                strain2_name = os.path.splitext(db2)[0]

                logger.debug(f'{"#"*100}\nGot {strain1_name} and {strain2_name}')
                if strain1_name >= strain2_name:
                    logger.debug(f'Skipping {strain1_name} >= {strain2_name}')
                    continue # no need to query strain against itself or a pair that was already seen
                logger.info(f'Querying {strain1_name} against {strain2_name}')

                query_aa_db = os.path.splitext(os.path.join(previous_pipeline_step_output_dir, db1))[0] + '_aa'

                target_aa_db = os.path.splitext(os.path.join(previous_pipeline_step_output_dir, db2))[0] + '_aa'

                aln_offsetted_db = os.path.join(pipeline_step_tmp_dir, f'{strain1_name}_vs_{strain2_name}.alnOffsettedDB')

                output_file_name = f'{strain1_name}_vs_{strain2_name}.m8'
                output_file_path = os.path.join(pipeline_step_output_dir, output_file_name)

                if num_of_aggregated_params > 0:
                    # params was already defined for this job batch. Save it before overridden
                    more_cmds.append(params)
                params = [query_aa_db,
                          target_aa_db,
                          aln_offsetted_db,
                          f'{pipeline_step_tmp_dir}/{strain1_name}_vs_{strain2_name}',
                          output_file_path]
                          #, f';!@#ls -1 {pipeline_step_output_dir} | grep "{strain1_name}_vs_{strain2_name}" | grep -v m8 | xargs rm']

                num_of_aggregated_params += 1
                if num_of_aggregated_params == num_of_cmds_per_job:
                    submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name,
                                         queue_name='pupkowebr -p -1', more_cmds=more_cmds, required_modules_as_list=[CONSTS.MMSEQS])
                    num_of_expected_results += 1
                    num_of_aggregated_params = 0
                    more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name,
                                 queue_name='pupkowebr -p -1', more_cmds=more_cmds, required_modules_as_list=[CONSTS.MMSEQS])
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=20)

    # send a subjob that removes (now instead of at the end) all mmseqs intermediate *FOLDERS* (e.g., 3465136234521948 etc..) in tmp_dir
    # does not remove (sge/cmds/log) files. only folders.
    # submit_pipeline_step(os.path.join(args.src_dir, 'remove_tmp_folders.py'), [previous_pipeline_step_output_dir], pipeline_step_tmp_dir,
    #                      job_name='remove_m8_files', queue_name=args.queue_name)

    # 4.	filter_blast.py
    # Input: (1) a path for a i_vs_j.tsv file (2) an output path (with a suffix as follows: i_vs_j_filtered.tsv. especially relevant for the wrapper).
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
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 100 if len(os.listdir(ORFs_dir)) > 100 else 50
    num_of_aggregated_params = 0
    more_cmds = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info('Filtering all vs all results...\n')
        logger.debug(f'Filtering the following files:\n' + '\n'.join(x for x in os.listdir(previous_pipeline_step_output_dir)))
        for blast_results_file in os.listdir(previous_pipeline_step_output_dir):
            fasta_file_prefix = os.path.splitext(blast_results_file)[0]
            output_file_name = f'{fasta_file_prefix}.{dir_name}'
            if num_of_aggregated_params > 0:  # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)

            params = [os.path.join(previous_pipeline_step_output_dir, blast_results_file),
                      os.path.join(pipeline_step_output_dir, output_file_name),
                      f'--identity_cutoff {args.identity_cutoff/100}', # needs to be normaized between 0 and 1
                      f'--e_value_cutoff {args.e_value_cutoff}']

            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name,
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params > 0:
            # don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name,
                                 queue_name=args.queue_name, more_cmds=more_cmds)
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=25)

    # send a subjob that removes all m8 files from the output dir as they are no longer needed and might cost a lot of inodes!
    # does not remove (sge/cmds/log) files. only folders.
    # submit_pipeline_step(os.path.join(args.src_dir, 'remove_tmp_folders.py'), [pipeline_step_tmp_dir], pipeline_step_tmp_dir,
    #                      job_name='remove_dirs_from_tmp', queue_name=args.queue_name)


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
        logger.info(f'Calling:\n{cmd}')
        subprocess.run(cmd, shell=True)
        # No need to wait...
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=30)


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
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing putative orthologs table...')
        job_name = os.path.split(script_path)[-1]
        params = [all_reciprocal_hits_file,
                  putative_orthologs_table_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name,
                             queue_name=args.queue_name)
        num_of_expected_results = 1
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=35)


    # 8   prepare_files_for_mcl.py
    # Input: (1) a path for a concatenated all reciprocal hits file (2) a path for a putative orthologs file (3) a path for an output folder
    # Output: an input file for MCL for each putative orthologs group
    # CANNOT be parallelized on cluster (if running on the concatenated file)
    step = '08'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_mcl_input_files'
    script_path = os.path.join(args.src_dir, 'prepare_files_for_mcl.py')
    num_of_expected_results = 0
    pipeline_step_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 25  # *times* the number of steps below. 250 in total per batch!
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs

    with open(os.path.join(os.path.split(putative_orthologs_table_path)[0], 'num_of_putative_sets.txt')) as f:
        num_of_putative_sets = int(f.read())

    if not os.path.exists(done_file_path):
        logger.info('Preparing files for MCL...')
        step = 10
        for i in range(1, num_of_putative_sets+1, step):
            if num_of_aggregated_params > 0:
                # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)
            first_mcl = str(i)
            last_mcl = str(min(i+step-1, num_of_putative_sets))  # 1-10, 11-20, etc... for step=10
            params = [all_reciprocal_hits_file,
                      putative_orthologs_table_path,
                      first_mcl,
                      last_mcl,
                      pipeline_step_output_dir]
            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=f'{dir_name}_{last_mcl}',  # {i+step-num_of_cmds_per_job*step}_
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=f'{dir_name}_{last_mcl}',
                                 queue_name=args.queue_name, more_cmds=more_cmds)
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=40)


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
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 100
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info('Executing MCL...')
        for putative_orthologs_group in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_group_prefix = os.path.splitext(putative_orthologs_group)[0]
            output_file_name = f'{putative_orthologs_group_prefix}.{dir_name}'
            if num_of_aggregated_params > 0:
                # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)
            params = [f'"{os.path.join(previous_pipeline_step_output_dir, putative_orthologs_group)}"',
                      f'"{os.path.join(pipeline_step_output_dir, output_file_name)}"']
            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name,
                                     queue_name=args.queue_name, required_modules_as_list=[CONSTS.MCL],
                                     more_cmds=more_cmds)
                num_of_expected_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=output_file_name,
                                 queue_name=args.queue_name, required_modules_as_list=[CONSTS.MCL],
                                 more_cmds=more_cmds)
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=45)


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
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 100
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info('Verifying clusters...')
        for putative_orthologs_group in os.listdir(previous_pipeline_step_output_dir):
            putative_orthologs_group_prefix = os.path.splitext(putative_orthologs_group)[0]
            job_name = os.path.split(putative_orthologs_group_prefix)[-1]
            if num_of_aggregated_params > 0:
                # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)
            params = [f'"{os.path.join(previous_pipeline_step_output_dir, putative_orthologs_group)}"',
                      f'"{os.path.join(pipeline_step_output_dir, putative_orthologs_group)}"']
            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name,
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name,
                                 queue_name=args.queue_name, more_cmds=more_cmds)
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=50)

    if num_of_expected_results == 0:
        error_msg = f'No ortholog groups were found in your dataset. Please try to lower the identity threshold (see Advanced options in the submission page) and re-submit your job.'
        fail(error_msg, error_file_path)

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
    phyletic_patterns_path = os.path.join(final_orthologs_table_path, 'phyletic_pattern.fas')
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Constructing final orthologs table...')
        job_name = os.path.split(final_orthologs_table_file_path)[-1]
        params = [putative_orthologs_table_path,
                  previous_pipeline_step_output_dir,
                  final_orthologs_table_file_path,
                  phyletic_patterns_path]
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=job_name, queue_name=args.queue_name)
        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)

        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=55)


    # 12.	extract_orthologs_sequences.py
    # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
    # Output: write the sequences of the orthologs group to the output file.
    step = '12'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_orthologs_groups_dna_sequences'
    script_path = os.path.join(args.src_dir, 'extract_orthologs_sequences.py')
    num_of_expected_results = 0
    orthologs_dna_sequences_dir_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    num_of_strains_path = os.path.join(final_orthologs_table_path, 'num_of_strains.txt')
    strains_names_path = os.path.join(final_orthologs_table_path, 'strains_names.txt')
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 250
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info(f'There are {len(os.listdir(previous_pipeline_step_output_dir))} verified clusters in {previous_pipeline_step_output_dir}.')
        logger.info('Extracting orthologs groups sequences according to final orthologs table...')
        logger.debug(f'The verified clusters in {previous_pipeline_step_output_dir} are the following:')
        logger.debug(os.listdir(previous_pipeline_step_output_dir))
        # og_number = 0
        with open(final_orthologs_table_file_path) as f:
            header_line = f.readline()
            first_delimiter_index = header_line.index(delimiter)
            final_table_header = header_line.rstrip()[first_delimiter_index + 1:] #remove "OG_name"
            for line in f:
                first_delimiter_index = line.index(delimiter)
                og_name = line[:first_delimiter_index]
                cluster_members = line.rstrip()[first_delimiter_index+1:] #remove "OG_name"
                output_file_name = og_name
                # og_number += 1
                if num_of_aggregated_params > 0:
                    # params was already defined for this job batch. Save it before overridden
                    more_cmds.append(params)
                params = [ORFs_dir,
                          f'"{final_table_header}"',  # should be flanked by quotes because it might contain spaces...
                          f'"{cluster_members}"',  # should be flanked by quotes because it might contain spaces...
                          os.path.join(orthologs_dna_sequences_dir_path, f'{output_file_name}_dna.fas')]
                num_of_aggregated_params += 1
                logger.debug(f'num_of_aggregated_params: {num_of_aggregated_params} of {num_of_cmds_per_job}')
                if num_of_aggregated_params == num_of_cmds_per_job:
                    submit_pipeline_step(script_path, params, pipeline_step_tmp_dir,
                                         job_name=output_file_name,
                                         queue_name=args.queue_name,
                                         more_cmds=more_cmds)
                    num_of_expected_results += 1
                    num_of_aggregated_params = 0
                    more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir,
                                 job_name=output_file_name,
                                 queue_name=args.queue_name,
                                 more_cmds=more_cmds)
            num_of_expected_results += 1

        # extract number of strains and their names for core genome analysis later on
        strains_names = final_table_header.split(delimiter)
        with open(strains_names_path, 'w') as f:
            f.write('\n'.join(strains_names)+'\n')

        num_of_strains = len(strains_names)
        with open(num_of_strains_path, 'w') as f:
            f.write(f'{num_of_strains}\n')

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
        # this step always needs to run!!
        # file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=60)


    # 13.  translate_fna_to_faa.py
    # Input: path to fna file and an faa file
    # Output: translate the fna to protein and write to the faa file
    # Can be parallelized on cluster
    step = '13'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_orthologs_groups_aa_sequences'
    script_path = os.path.join(args.src_dir, 'translate_fna_to_faa.py')
    num_of_expected_results = 0
    orthologs_aa_sequences_dir_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 250
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info('Translating orthologs groups sequences...')
        for fasta_file in os.listdir(orthologs_dna_sequences_dir_path):
            file_path = os.path.join(orthologs_dna_sequences_dir_path, fasta_file)
            output_path = os.path.join(orthologs_aa_sequences_dir_path, fasta_file.replace('_dna.fas', '_aa.fas'))

            job_name = f'{os.path.splitext(fasta_file)[0]}_to_aa'
            if num_of_aggregated_params > 0:  # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)
            params = [file_path,
                      output_path]
            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name,
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name,
                                 queue_name=args.queue_name, more_cmds=more_cmds)
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)

        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=65)


    # 14.	align_orthologs_group.py
    # Input: (1) A path to an unaligned amino acid sequences file (2) An output file path
    # Output: aligned sequences
    # Can be parallelized on cluster
    step = '14'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_aligned_aa_orthologs_groups'
    script_path = os.path.join(args.src_dir, 'align_orthologs_group.py')
    num_of_expected_results = 0
    aa_alignments_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 250
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    if not os.path.exists(done_file_path):
        logger.info('Aligning orthologs groups...')
        for og_file in os.listdir(orthologs_aa_sequences_dir_path):
            # if not og_file.endswith('fas'):
            #     # why should that happen?
            #     continue
            og_path = os.path.join(orthologs_aa_sequences_dir_path, og_file)
            og_file_prefix = os.path.splitext(og_file)[0]
            alignment_path = os.path.join(aa_alignments_path, f'{og_file_prefix}_mafft.fas')
            if num_of_aggregated_params > 0:
                # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)
            params = [og_path, alignment_path]
            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=og_file_prefix,
                                     queue_name=args.queue_name, required_modules_as_list=[CONSTS.MAFFT],
                                     more_cmds=more_cmds)
                num_of_expected_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=og_file_prefix,
                                 queue_name=args.queue_name, required_modules_as_list=[CONSTS.MAFFT],
                                 more_cmds=more_cmds)
            num_of_expected_results += 1

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results,
                         error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=70)


    #15.	extract aligned core genome.py IN FACT, THE CORE PROTEOME IS EXTRACTED
    step = '15'
    logger.info(f'Step {step}: {"_" * 100}')
    dir_name = f'{step}_aligned_core_proteome'
    script_path = os.path.join(args.src_dir, 'extract_core_genome.py')
    num_of_expected_results = 1
    aligned_core_proteome_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    aligned_core_proteome_file_path = os.path.join(aligned_core_proteome_path, 'aligned_core_proteome.fas')
    core_ogs_names_file_path = os.path.join(aligned_core_proteome_path, 'core_ortholog_groups_names.txt')
    core_length_file_path = os.path.join(aligned_core_proteome_path, 'core_length.txt')
    with open(num_of_strains_path) as f:
        num_of_strains = f.read().rstrip()
    if not os.path.exists(done_file_path):
        logger.info('Extracting aligned core proteome...')

        params = [aa_alignments_path, num_of_strains, strains_names_path,
                  aligned_core_proteome_file_path,
                  core_ogs_names_file_path,
                  core_length_file_path,
                  f'--core_minimal_percentage {args.core_minimal_percentage}']  # how many members induce a core group?
        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name='core_proteome',
                             queue_name=args.queue_name)

        wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results,
                         error_file_path)
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=75)


    # 16.	reconstruct_species_phylogeny.py
    step = '16'
    logger.info(f'Step {step}: {"_" * 100}')
    dir_name = f'{step}_species_phylogeny'
    script_path = os.path.join(args.src_dir, 'reconstruct_species_phylogeny.py')
    num_of_expected_results = 1
    phylogeny_path, phylogeny_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    phylogenetic_raw_tree_path = os.path.join(phylogeny_path, 'final_species_tree.txt')
    # no need to wait now. Wait before plotting the tree!
    # done_file_path = os.path.join(done_files_dir, f'{step}_reconstruct_species_phylogeny.txt')
    start_tree = time()
    if not os.path.exists(phylogenetic_raw_tree_path):
        logger.info('Reconstructing species phylogeny...')

        if 'pupkoweb' not in args.queue_name:
            num_of_cpus = 20  # anyhow it can't be more than 20! o.w., "qsub: Job violates queue and/or server resource limits"
        else:
            # don't use more than 5 cpus when running on pupkoweb
            num_of_cpus = 5
        params = [aligned_core_proteome_file_path,
                  phylogenetic_raw_tree_path,
                  f'--cpu {num_of_cpus}']  # both the Q and RAxML should get this as a parameter
                  # '--model PROTGAMMAILG',
                  # f'--num_of_bootstrap_iterations {100 if args.bootstrap=="yes" else 0}',
        if args.outgroup:
            if args.outgroup in strains_names:
                params += [f'--outgroup {args.outgroup}']
            else:
                logger.info(f'Outgroup {args.outgroup} was specified but it is not one of the input species:\n'
                            f'{",".join(sorted(strains_names))}\nAn unrooted tree is going to be reconstructed')
        if args.bootstrap == 'yes' and number_of_genomes < 120:  # allow bootstrap only for less than 120 genomes
            params += ['--num_of_bootstrap_iterations 100']

        submit_pipeline_step(script_path, params, phylogeny_tmp_dir, job_name='tree_reconstruction',
                             queue_name=args.queue_name, required_modules_as_list=[CONSTS.RAXML], num_of_cpus=num_of_cpus)
        # no need to wait now. Wait before plotting the tree!
        # Still, in order to get more than one thread, allow few seconds to the subjobs to be submitted
        # (before the next batch is submitted and pounds the cluster)
        sleep(60)
    else:
        logger.info(f'Raw tree file {phylogenetic_raw_tree_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=80)


    # 17.	extract_orfs_statistics.py
    # Input: (1) A path to ORFs file (2) An output path to ORFs counts (3) An output path to GC content
    # Output: write the number of ORFs and GC content to the output files (respectively)
    # Can be parallelized on cluster
    step = '17'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_orfs_statistics'
    script_path = os.path.join(args.src_dir, 'extract_orfs_statistics.py')
    num_of_expected_orfs_results = 0
    orfs_statistics_path, orfs_statistics_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 100
    num_of_aggregated_params = 0
    more_cmds = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    start_orf_stats = time()
    if not os.path.exists(done_file_path):
        logger.info('Collecting orfs counts...')

        for file in os.listdir(ORFs_dir):
            orf_path = os.path.join(ORFs_dir, file)
            strain_name = os.path.splitext(file)[0]
            orfs_count_output_path = os.path.join(orfs_statistics_path, f'{strain_name}.orfs_count')
            gc_content_output_path = os.path.join(orfs_statistics_path, f'{strain_name}.gc_content')

            if num_of_aggregated_params > 0:
                # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)
            # args.orf_path, args.orfs_counts_output_path, args.orfs_gc_output_path
            params = [orf_path, orfs_count_output_path, gc_content_output_path]
            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, orfs_statistics_tmp_dir,
                                     job_name=f'{strain_name}_orfs_stats',
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_orfs_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params > 0:
            # don't forget the last batch!!
            submit_pipeline_step(script_path, params, orfs_statistics_tmp_dir,
                                 job_name=f'{strain_name}_orfs_stats',
                                 queue_name=args.queue_name, more_cmds=more_cmds)
            num_of_expected_orfs_results += 1

        # no need to wait now. Wait before plotting the statistics!

        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    #edit_progress(output_html_path, progress=85)


    # 18.	induce_dna_msa_by_aa_msa.py
    # Input: (1) An aa alignment (2) An unaligned dna file (3) An output file path
    # Output: write the codon alignment induced by the aa alignment to the output file
    # Can be parallelized on cluster
    step = '18'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_induce_dna_msa_by_aa_msa'
    script_path = os.path.join(args.src_dir, 'induce_dna_msa_by_aa_msa.py')
    num_of_expected_induced_results = 0
    dna_alignments_path, induced_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    num_of_cmds_per_job = 250
    num_of_aggregated_params = 0
    more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    start_induced = time()
    if not os.path.exists(done_file_path):
        logger.info(f'Inducing dna alignments...\n(from {aa_alignments_path})')

        for og_file in os.listdir(aa_alignments_path):
            aa_alignment_path = os.path.join(aa_alignments_path, og_file)
            dna_unaligned_path = os.path.join(orthologs_dna_sequences_dir_path, og_file.replace('aa_mafft', 'dna'))
            dna_induced_alignment_path = os.path.join(dna_alignments_path, og_file.replace('aa_mafft','dna_induced'))
            if num_of_aggregated_params > 0:
                # params was already defined for this job batch. Save it before overridden
                more_cmds.append(params)
            params = [aa_alignment_path, dna_unaligned_path, dna_induced_alignment_path]
            num_of_aggregated_params += 1
            if num_of_aggregated_params == num_of_cmds_per_job:
                submit_pipeline_step(script_path, params, induced_tmp_dir, job_name=f'induced_{og_file}',
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_induced_results += 1
                num_of_aggregated_params = 0
                more_cmds = []

        if num_of_aggregated_params>0:
            #don't forget the last batch!!
            submit_pipeline_step(script_path, params, induced_tmp_dir, job_name=f'induced_{og_file}',
                                 queue_name=args.queue_name, more_cmds=more_cmds)
            num_of_expected_induced_results += 1

        # no need to wait now. Wait before moving the results dir!
        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 19.	extract_groups_sizes_frequency
    step = '19'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_groups_sizes_frequency'
    group_sizes_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    groups_sizes_frequency_file_prefix = os.path.join(group_sizes_path, 'groups_sizes_frequency')
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Collecting sizes...')

        group_sizes = []
        with open(final_orthologs_table_file_path) as f:
            final_table_header = f.readline().rstrip()
            for line in f:
                cluster_members = line.rstrip()
                size = sum(bool(item) for item in cluster_members.split(delimiter))  # count non empty entries
                size -= 1 # don't count group name
                group_sizes.append(str(size))

        groups_sizes_frequency_raw_file_path = groups_sizes_frequency_file_prefix + '.txt'
        with open(groups_sizes_frequency_raw_file_path, 'w') as f:
            f.write('\n'.join(group_sizes)) #f.write('\n'.join([f'{size},{group_size_to_counts_dict[size]}' for size in group_size_to_counts_dict]))

        groups_sizes_frequency_png_file_path = groups_sizes_frequency_file_prefix + '.png'
        generate_bar_plot(groups_sizes_frequency_raw_file_path, groups_sizes_frequency_png_file_path,
            xlabel='Orthologous group size', ylabel='Count')

        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


    # 20.	plot_orfs_statistics
    step = '20'
    logger.info(f'Step {step}: {"_"*100}')
    dir_name = f'{step}_orfs_plots'
    orfs_plots_path, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
    orfs_counts_frequency_file = os.path.join(orfs_plots_path, 'orfs_counts.txt')
    orfs_gc_content_file = os.path.join(orfs_plots_path, 'orfs_gc_contents.txt')
    done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
    if not os.path.exists(done_file_path):

        wait_for_results('extract_orfs_statistics.py', orfs_statistics_tmp_dir,
                         num_of_expected_orfs_results, error_file_path=error_file_path, start=start_orf_stats)

        logger.info('Concatenating orfs counts...')
        cmd = f'cat {orfs_statistics_path}/*.orfs_count > {orfs_counts_frequency_file}'
        subprocess.run(cmd, shell = True)

        logger.info('Concatenating orfs gc contents...')
        cmd = f'cat {orfs_statistics_path}/*.gc_content > {orfs_gc_content_file}'
        subprocess.run(cmd, shell = True)

        # No need to wait...

        logger.info('Ploting violines...')
        orfs_counts_frequency_png_file_path = orfs_counts_frequency_file.replace('txt', 'png')
        generate_boxplot(orfs_counts_frequency_file, orfs_counts_frequency_png_file_path,
                         xlabel='\nORF count per genome', dpi=300)

        orfs_gc_content_png_file_path = orfs_gc_content_file.replace('txt', 'png')
        generate_boxplot(orfs_gc_content_file, orfs_gc_content_png_file_path,
                         xlabel='\nGC content per genome', dpi=300)

        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')
    edit_progress(output_html_path, progress=85)


    # # 21.	plot species phylogeny
    # step = '21'
    # logger.info(f'Step {step}: {"_" * 100}')
    # dir_name = f'{step}_plot_tree'
    # done_file_path = os.path.join(done_files_dir, f'{step}_plot_species_phylogeny.txt')
    # if int(num_of_strains) > 3 and not os.path.exists(done_file_path):
    #     logger.info('Ploting species phylogeny...')
    #
    #     # wait for the raw tree here
    #     wait_for_results('reconstruct_species_phylogeny.py', phylogeny_tmp_dir,
    #                      num_of_expected_results=1, error_file_path=error_file_path)
    #
    #     phylogenetic_png_tree_path = phylogenetic_raw_tree_path.replace('txt', 'png')
    #     try:
    #         generate_tree_plot(phylogenetic_raw_tree_path, phylogenetic_png_tree_path)
    #     except:
    #         if not os.path.exists(phylogenetic_raw_tree_path):
    #             logger.fatal(f'Raw tree file {phylogenetic_raw_tree_path} does not exist! Possible bug?')
    #         else:
    #             logger.fatal(f'Raw tree file {phylogenetic_raw_tree_path} exists but could not generate a png! Check RAxML logs for further info.')
    #
    #     file_writer.write_to_file(done_file_path, '.')
    # else:
    #     logger.info(f'Number of strains is {num_of_strains}')
    #     logger.info(f'Number of strains < 4 or done file {done_file_path} already exists.\nSkipping step...')


    # Final step: gather relevant results, zip them together and update html file
    logger.info(f'FINAL STEP: {"_"*100}')
    final_output_dir_name = f'{CONSTS.WEBSERVER_NAME}_{run_number}_outputs'
    final_output_dir, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, final_output_dir_name)
    done_file_path = os.path.join(done_files_dir, f'{final_output_dir_name}.txt')
    if not os.path.exists(done_file_path):
        logger.info('Gathering results to final output dir...')

        # move ORFs folder
        add_results_to_final_dir(ORFs_dir, final_output_dir, copy=True)

        # move orthologs table
        add_results_to_final_dir(final_orthologs_table_path, final_output_dir, copy=True)

        # move unaligned dna sequences
        add_results_to_final_dir(orthologs_dna_sequences_dir_path, final_output_dir)

        # move unaligned aa sequences
        add_results_to_final_dir(orthologs_aa_sequences_dir_path, final_output_dir)

        # move aligned aa sequences
        add_results_to_final_dir(aa_alignments_path, final_output_dir)

        # move groups sizes
        add_results_to_final_dir(group_sizes_path, final_output_dir)

        # move orfs statistics
        add_results_to_final_dir(orfs_statistics_path, final_output_dir)

        # move orfs plot
        add_results_to_final_dir(orfs_plots_path, final_output_dir)

        wait_for_results('induce_dna_msa_by_aa_msa.py', induced_tmp_dir,
                         num_of_expected_results=num_of_expected_induced_results,
                         error_file_path=error_file_path, start=start_induced)
        # move induced dna sequences
        add_results_to_final_dir(dna_alignments_path, final_output_dir)

        # wait for the raw tree here
        wait_for_results('reconstruct_species_phylogeny.py', phylogeny_tmp_dir,
                         num_of_expected_results=1, error_file_path=error_file_path, start=start_tree)
        edit_progress(output_html_path, progress=98)

        # move core proteome dir
        add_results_to_final_dir(aligned_core_proteome_path, final_output_dir)

        # move species tree dir
        add_results_to_final_dir(phylogeny_path, final_output_dir)

        logger.info('Zipping results folder...')
        shutil.make_archive(final_output_dir, 'zip', final_output_dir)

        logger.info(f'Moving results to parent dir... ({meta_output_dir})')
        try:
            shutil.move(f'{final_output_dir}.zip', meta_output_dir)
        except shutil.Error as e:
            logger.error(e.args[0])
        try:
            shutil.move(final_output_dir, meta_output_dir)
        except shutil.Error as e:
            logger.error(e.args[0])

        file_writer.write_to_file(done_file_path, '.')
    else:
        logger.info(f'done file {done_file_path} already exists.\nSkipping step...')

    logger.info('Editing results html...')
    edit_success_html(output_html_path, meta_output_dir, final_output_dir_name, run_number, CONSTS)

    edit_progress(output_html_path, progress=100, active=False)

    status = 'is done'

except Exception as e:
    status = 'was failed'
    from time import ctime
    import os
    import logging
    logger = logging.getLogger('main')  # use logger instead of printing

    error_msg = f'{CONSTS.WEBSERVER_NAME} failed :('
    if os.path.exists(error_file_path):
        with open(error_file_path) as error_f:
            error_txt = error_f.read()
            logger.error(f'error.txt file says:\n{error_txt}')
            error_msg = f'The job was failed due to the following reason:<br>{error_txt}'

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    logger.error(f'\n\n{"$" * 100}\n\n{error_msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\ne: {e}\n\n{"$" * 100}')

    edit_failure_html(output_html_path, run_number, error_msg, CONSTS)
    edit_progress(output_html_path, active=False)

    notify_admin(meta_output_dir, meta_output_url, run_number, CONSTS)

end = time()

results_location = output_url if remote_run else args.output_dir
msg = f'M1CR0B1AL1Z3R pipeline {status}'
if status == 'is done':
    msg += f' (Took {measure_time(int(end-start))}).\nResults can be found at {results_location}.\nPlease note that the results will be kept in the server for three months.'
else:
    msg += f'. For further information please visit: {results_location}'
logger.info(msg)

logger.info(f'Sending a notification email to {args.email}')
try:
    send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>', args.email, subject=f'{CONSTS.WEBSERVER_NAME} run number {run_number} {status}.', content=msg)
except:
    logger.info(f'\nFailed sending notification to {args.email}\n')


logger.info('Cleaning up...')
if status == 'is done':

    if 'oren' in args.email:
        # 21.	extract_promoters_and_orfs
        # Input: (1) A path to a genome (2) Prodigal output file with ORFs coordinates
        # Output: A fasta file containing promoters+genes (all ORFS of the given coordinates + k[=300] upstream bases)
        # Can be parallelized on cluster
        step = '21'
        logger.info(f'Step {step}: {"_"*100}')
        dir_name = f'{step}_extract_promoters_and_orfs'
        script_path = f'{args.src_dir}/../post_analysis/extract_promoters_and_orfs.py'
        num_of_expected_results = 0
        pipeline_step_output_dir_21, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
        done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
        num_of_cmds_per_job = 5
        num_of_aggregated_params = 0
        more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        if not os.path.exists(done_file_path):
            logger.info('Extracting promoters...')
            for genome_file in os.listdir(data_path):
                genome_path = os.path.join(data_path, genome_file)
                genome_file_prefix = os.path.splitext(genome_file)[0]
                orfs_path = os.path.join(ORFs_dir, f'{genome_file_prefix}.01_ORFs')
                output_file_name = f'{genome_file_prefix}.promoter_and_orf'
                if num_of_aggregated_params > 0:  # params was already defined for this job batch. Save it before overridden
                    more_cmds.append(params)

                output_prefix = os.path.join(pipeline_step_output_dir_21, output_file_name)
                params = [genome_path,
                          orfs_path,
                          output_prefix,
                          f'--promoters_length {300}'] # (optional)

                num_of_aggregated_params += 1
                if num_of_aggregated_params == num_of_cmds_per_job:
                    submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=genome_file_prefix,
                                         queue_name=args.queue_name, more_cmds=more_cmds)
                    num_of_expected_results += 1
                    num_of_aggregated_params = 0
                    more_cmds = []

            if num_of_aggregated_params>0:
                #don't forget the last batch!!
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=genome_file_prefix,
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_results += 1

            wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
            file_writer.write_to_file(done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists.\nSkipping step...')

        # 22.	extract_orthologs_sequences.py
        # Input: (1) a row from the final orthologs table (2) a path for a directory where the genes files are at (3) a path for an output file.
        # Output: write the sequences of the orthologs group to the output file.
        # Can be parallelized on cluster
        step = '22'
        logger.info(f'Step {step}: {"_"*100}')
        dir_name = f'{step}_orthologs_groups_dna_sequences'
        script_path = os.path.join(args.src_dir, 'extract_orthologs_sequences.py')
        num_of_expected_results = 0
        pipeline_step_output_dir_22, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
        done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
        num_of_cmds_per_job = 100
        num_of_aggregated_params = 0
        more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        if not os.path.exists(done_file_path):
            logger.info(f'There are {len(os.listdir(previous_pipeline_step_output_dir))} verified clusters in {previous_pipeline_step_output_dir}.')
            logger.info('Extracting orthologs groups sequences according to final orthologs table...')
            logger.debug(f'The verified clusters in {previous_pipeline_step_output_dir} are the following:')
            logger.debug(os.listdir(pipeline_step_output_dir_21))
            # og_number = 0
            with open(final_orthologs_table_file_path) as f:
                header_line = f.readline()
                first_delimiter_index = header_line.index(delimiter)
                final_table_header = header_line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"
                for line in f:
                    first_delimiter_index = line.index(delimiter)
                    og_name = line[:first_delimiter_index]
                    cluster_members = line.rstrip()[first_delimiter_index + 1:]  # remove "OG_name"
                    output_file_name = og_name
                    # og_number += 1
                    if num_of_aggregated_params > 0:
                        # params was already defined for this job batch. Save it before overridden
                        more_cmds.append(params)
                    params = [pipeline_step_output_dir_21,
                              f'"{final_table_header}"',  # should be flanked by quotes because it might contain spaces...
                              f'"{cluster_members}"',  # should be flanked by quotes because it might contain spaces...
                              os.path.join(pipeline_step_output_dir_22, f'{output_file_name}_dna.fas')]
                    num_of_aggregated_params += 1
                    logger.debug(f'num_of_aggregated_params: {num_of_aggregated_params} of {num_of_cmds_per_job}')
                    if num_of_aggregated_params == num_of_cmds_per_job:
                        submit_pipeline_step(script_path, params, pipeline_step_tmp_dir,
                                             job_name=output_file_name,
                                             queue_name=args.queue_name,
                                             more_cmds=more_cmds)
                        num_of_expected_results += 1
                        num_of_aggregated_params = 0
                        more_cmds = []

            if num_of_aggregated_params>0:
                #don't forget the last batch!!
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir,
                                     job_name=output_file_name,
                                     queue_name=args.queue_name,
                                     more_cmds=more_cmds)
                num_of_expected_results += 1

            wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results, error_file_path)
            file_writer.write_to_file(done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists.\nSkipping step...')


        # 23.	align_orthologs_group.py
        # Input: (1) A path to an unaligned sequences file (2) An output file path
        # Output: aligned sequences
        step = '23'
        logger.info(f'Step {step}: {"_"*100}')
        dir_name = f'{step}_aligned_dna_orthologs_groups_with_promoter'
        script_path = os.path.join(args.src_dir, 'align_orthologs_group.py')
        num_of_expected_results = 0
        pipeline_step_output_dir_23, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
        done_file_path = os.path.join(done_files_dir, f'{dir_name}.txt')
        num_of_cmds_per_job = 100
        num_of_aggregated_params = 0
        more_cmds = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        if not os.path.exists(done_file_path):
            logger.info('Aligning orthologs groups...')
            for og_file in os.listdir(pipeline_step_output_dir_22):
                og_path = os.path.join(pipeline_step_output_dir_22, og_file)
                og_file_prefix = os.path.splitext(og_file)[0]
                alignment_path = os.path.join(pipeline_step_output_dir_23, f'{og_file_prefix}_mafft.fas')
                if num_of_aggregated_params > 0:
                    # params was already defined for this job batch. Save it before overridden
                    more_cmds.append(params)
                params = [og_path, alignment_path, '--type nuc']
                num_of_aggregated_params += 1
                if num_of_aggregated_params == num_of_cmds_per_job:
                    submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=og_file_prefix,
                                         queue_name=args.queue_name, required_modules_as_list=[CONSTS.MAFFT],
                                         more_cmds=more_cmds)
                    num_of_expected_results += 1
                    num_of_aggregated_params = 0
                    more_cmds = []

            if num_of_aggregated_params>0:
                #don't forget the last batch!!
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=og_file_prefix,
                                     queue_name=args.queue_name, required_modules_as_list=[CONSTS.MAFFT],
                                     more_cmds=more_cmds)
                num_of_expected_results += 1

            wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results,
                             error_file_path)
            file_writer.write_to_file(done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists.\nSkipping step...')

        # 24.	adjust_tree_to_msa.py
        # Input: (1) A path to an MSA file to which the tree should be adjusted
        #        (2) A path to a species tree that contains (at least) all the species in the input MSA
        #        (3) A path to a folder in which a txt file (with the same name as the msa_file) will be created. All the species that do not appear in the msa (and thus will be removed) will be written to the file that was created',
        #        (4) A path to a file in which the prunned tree will be written
        # Output: an adjusted tree per msa at (4)
        step = '24'
        logger.info(f'Step {step}: {"_"*100}')
        dir_name = f'{step}_prunned_trees'
        script_path = f'{args.src_dir}/../post_analysis/adjust_tree_to_msa.py'
        num_of_expected_results = 0
        pipeline_step_output_dir_24, pipeline_step_tmp_dir = prepare_directories(args.output_dir, tmp_dir, dir_name)
        done_file_path = os.path.join(done_files_dir, f'{step}_{dir_name}.txt')
        num_of_cmds_per_job = 250
        num_of_aggregated_params = 0
        more_cmds = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        if not os.path.exists(done_file_path):
            phylogenetic_raw_tree_path = phylogenetic_raw_tree_path.replace('outputs', final_output_dir_name)  # the tree was moved to the final dir..
            phylogenetic_raw_tree_path_without_bootstrap_values = phylogenetic_raw_tree_path.replace('.txt', '_no_bootstrap.txt')
            logger.info('Removing bootstrap values from tree...')
            remove_bootstrap_values(phylogenetic_raw_tree_path, phylogenetic_raw_tree_path_without_bootstrap_values)

            logger.info('Prunning trees...')
            for msa_file in os.listdir(pipeline_step_output_dir_23):
                msa_path = os.path.join(pipeline_step_output_dir_23, msa_file)
                msa_file_prefix = os.path.splitext(msa_file)[0]
                prunned_tree_path = os.path.join(pipeline_step_output_dir_24, f'{msa_file_prefix}.tree')
                if num_of_aggregated_params > 0:
                    # params was already defined for this job batch. Save it before overridden
                    more_cmds.append(params)
                params = [msa_path, phylogenetic_raw_tree_path_without_bootstrap_values, pipeline_step_tmp_dir, prunned_tree_path]
                num_of_aggregated_params += 1
                if num_of_aggregated_params == num_of_cmds_per_job:
                    submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=msa_file_prefix,
                                         queue_name=args.queue_name, more_cmds=more_cmds)
                    num_of_expected_results += 1
                    num_of_aggregated_params = 0
                    more_cmds = []

            if num_of_aggregated_params > 0:
                # don't forget the last batch!!
                submit_pipeline_step(script_path, params, pipeline_step_tmp_dir, job_name=msa_file_prefix,
                                     queue_name=args.queue_name, more_cmds=more_cmds)
                num_of_expected_results += 1

            wait_for_results(os.path.split(script_path)[-1], pipeline_step_tmp_dir, num_of_expected_results,
                             error_file_path)
            file_writer.write_to_file(done_file_path, '.')
        else:
            logger.info(f'done file {done_file_path} already exists.\nSkipping step...')

        # preventing folder's deletion (output dir is being deleted
        logger.info(f'Moving {pipeline_step_output_dir_24} TO {os.path.join(meta_output_dir, dir_name)}')
        try:
            shutil.move(pipeline_step_output_dir_23, meta_output_dir)
            shutil.move(pipeline_step_output_dir_24, meta_output_dir)
        except FileExistsError:
            pass

    if remote_run and run_number.lower() != 'example' and 'oren' not in args.email:

        # remove intermediate results (including tmp_dir)
        remove_path(args.output_dir)

        # remove raw data from the server
        for path_to_remove in [os.path.join(meta_output_dir, final_output_dir_name, x) for x in
                               ['12_orthologs_groups_dna_sequences',
                                '13_orthologs_groups_aa_sequences',
                                '14_aligned_aa_orthologs_groups',
                                '15_aligned_core_proteome',
                                '17_orfs_statistics',
                                '18_induce_dna_msa_by_aa_msa']]:

            # remove intermediate results
            remove_path(path_to_remove)

# remove data
try:
    remove_path(unzipped_data_path)
except:
    pass
logger.info('Done.')
