#!/data/shared/python/anaconda3-5.1.0/bin/python3.6

import os
import sys
import sh
import cgi
import cgitb
import subprocess
from time import time, ctime
from random import randint

# sys.path.append('/bioseq/bioSequence_scripts_and_constants')
sys.path.append('/bioseq/microbializer/auxiliaries')
import MICROBIALIZER_CONSTANTS as CONSTS
from directory_creator import create_dir
from email_sender import send_email


def write_to_debug_file(cgi_debug_path, content):
    with open(cgi_debug_path, 'a') as f:
        f.write(f'{ctime()}: {content}\n')


def write_html_prefix(output_path, run_number):
    with open(output_path, 'w') as f:
        f.write(f'''<html><head>

    <meta http-equiv="cache-control" content="no-cache, must-revalidate, post-check=0, pre-check=0" />
    <meta http-equiv="cache-control" content="max-age=0" />
    <meta http-equiv="expires" content="0" />
    <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
    <meta http-equiv="pragma" content="no-cache" />
    {CONSTS.RELOAD_TAGS}

    <title>MICROBIALIZER Job {run_number}</title>
    <link rel="shortcut icon" type="image/x-icon" href="{CONSTS.MICROBIALIZER_URL}/pics/logo.gif" />

    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
    <link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">

    <link rel="stylesheet" href="{CONSTS.MICROBIALIZER_URL}/css/general.css">

    </head><body>
    <nav role="navigation" class="navbar navbar-fixed-top">
        <div class="jumbotron" id="jumbo">
            <div class="container">            
                <div class="row" id="title-row" align="center">
                    <div class="col-md-1">
                    </div>
                    <div class="col-md-10">
                        <span id="server-title">M1CR0B1AL1Z3R</span>
                        <img src="{CONSTS.MICROBIALIZER_URL}/pics/logo.gif" id="nav_bar_image" class="img-rounded">
                        <br><span id="sub-title">Analyze bacterial genomic sequences. Easily.</span>
                    </div>
                </div>
            </div>       
        </div>
    </nav>
    <div id="behind-nav-bar-results">
    </div>
''')


def write_info_paragraph_to_html(output_path):
    with open(output_path, 'a') as f:
        f.write(f"""<br><div class="container" style="font-size: 20px;" align="justify"> 
<H1 align=center>Job Status - <FONT color='red'>RUNNING</FONT></h1>
<br>MICROBIALIZER is now processing your request. This page will be automatically updated every {CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. A link to this page was sent to your email in case you wish to view these results at a later time without recalculating them. Please note that the results will be kept in the server for three months.
</div>
<br>""")


def write_running_parameters_to_html(output_path, job_title, number_of_duplicates, urls_to_reads_files, files_names,
                                     MMU, len_threshold, qlty_threshold, number_of_clones_to_analyze, raw_data_suffix,
                                     mass_spec_seq, lib_file_name, f_umi, r_umi):
    with open(output_path, 'a') as f:

        # regular params row
        f.write("""<div class="container"><u><h3>Running Parameters:</h3></u><br>""")

        f.write('<div class="row" style="font-size: 20px;">')
        if job_title != '':
            f.write('<div class="col-md-6">')
            f.write(f'<b>Job title: </b>{job_title}<br><br>')
            f.write('</div>')
            f.write('</div><div class="row" style="font-size: 20px;"><br>')


        # show only path suffix... http://MICROBIALIZER.tau.ac.il/results/1520442230/reads/run1/242_R1.fastq
        f.write(f'''<div class="col-md-4">
                        <b>Fasta folder:</b><br>
                        <a href="{urls_to_reads_files[i * 2]}" target=_blank>{files_names[i*2]}</A><br>
                    </div>''')
        f.write('</div>')

        f.write(f'''
    <br><u><h3>Advanced Parameters:</h3></u><br>
    <div class="row">
        <div class="col-md-2">
            <b>Min read length: </b>{len_threshold}<br>
        </div>
        <div class="col-md-2">
            <b>Min read quality: </b>{qlty_threshold}<br>
        </div>
        <div class="col-md-2">
            <b>Raw data format: </b>{raw_data_suffix}<br>
        </div>
        <div class="col-md-6">
            <b>reference library: </b>{lib_file_name}<br>
        </div>
    </div>

</div>''')


# def write_pair_file(debug_path, pair, run_content, run_filename, run_dir):
#     if run_filename.endswith('gz'):
#         local_file_name = f'R{pair}.fastq.gz'
#     else:
#         local_file_name = f'R{pair}.fastq'
#     open_operator = open
#     with open(debug_path, 'a') as f:
#         f.write(f'{run_filename} ({local_file_name}) is being handled with {open_operator}\n')
#     local_file_path = os.path.join(run_dir, local_file_name)
#     with open_operator(local_file_path, 'wb') as f:
#         f.write(run_content)
#
#     # avoid double zipping:
#     if local_file_path.endswith('gz'):
#         sh.gunzip(local_file_path)
#
#     with open(debug_path, 'a') as f:
#         f.write(f'R{pair} was handled successfully\n')


# def process_uploaded_file(run_dir, debug_path):
#     # for debugging
#     with open(debug_path, 'a') as f:
#         f.write(f'{"#"*80}\nUploading file of {run}\n')
#     run_R1_filename = form[f'{run}_R1'].filename
#     run_R2_filename = form[f'{run}_R2'].filename
#     # for debugging
#     with open(debug_path, 'a') as f:
#         f.write(
#             f'file names are:\n{run_R1_filename} (of type {type(form[run + "_R1"].value)}) and {run_R2_filename} (of type {type(form[run + "_R2"].value)})\n')
#     run_R1_content = form[f'{run}_R1'].value
#     run_R2_content = form[f'{run}_R2'].value
#     # for debugging
#     with open(debug_path, 'a') as f:
#         f.write(
#             f'{run_R1_filename} first 100 chars are: {run_R1_content[:100]}\n{run_R2_filename} first 100 chars are: {run_R2_content[:100]}\n')
#
#     write_pair_file(debug_path, '1', run_R1_content, run_R1_filename, run_dir)
#
#     write_pair_file(debug_path, '2', run_R2_content, run_R2_filename, run_dir)
#
#     return run_R1_filename, run_R2_filename


def prepare_parameters_file(parameters_file, wd):
    with open(parameters_file, 'w') as f:
        f.write(f"""
    #String that represents the path to a folder where the output dir will be generated
    {wd}
    """)


def write_cmds_file(cmds_file_path, data_dir, output_dir, user_email, run_number):
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
    new_line_delimiter = ';!@#'
    # the code contains features that are exclusive to Python3.6 (or higher)!
    # repseqio and mixcr require java 1.8 (or higher)
    required_modules = ' '.join(
        ['python/anaconda_python-3.6.4', 'mafft/mafft7313', 'MCL-edge/mcl-14-137'])
    with open(cmds_file_path, 'w') as f:
        f.write(f'module load {required_modules}')
        f.write(new_line_delimiter)
        f.write(f'{" ".join(["python", CONSTS.MAIN_SCRIPT, data_dir, output_dir, user_email])}\tMICROBIALIZER_{run_number}')


# prints detailed error report on BROWSER when cgi crashes
# This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
cgitb.enable()

# print_hello_world() # for debugging
form = cgi.FieldStorage()  # extract POSTed object

# random_chars = "".join(choice(string.ascii_letters + string.digits) for x in range(20))
run_number = str(round(time())) + str(
    randint(10 ** 19, 10 ** 20 - 1))  # adding 20 random digits to prevent users see data that are not their's
if True:
    run_number = 'debug'  # str(round(time())) + str(randint(1000,9999)) # adding 4 random figures to prevent users see data that are not their's

results_url = os.path.join(CONSTS.MICROBIALIZER_RESULTS_URL, run_number)
output_url = os.path.join(results_url, 'output.html')

wd = os.path.join(CONSTS.MICROBIALIZER_RESULTS_DIR, run_number)
create_dir(wd)
output_html_path = os.path.join(wd, 'output.html')
cgi_debug_path = os.path.join(wd, 'cgi_debug.txt')

# print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
# print_hello_world() # comment out for debugging
# print_hello_world(output_html_path, run_number) # comment out for debugging

write_html_prefix(output_html_path, run_number)  # html's prefix must be written BEFORE redirecting...

print(f'Location: {output_url}')  # Redirects to the results url. MUST appear before any other print.
print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
sys.stdout.flush()  # must be flushed immediately!!!

# Send me a notification email every time there's a new request
send_email(smtp_server=CONSTS.SMTP_SERVER, sender=CONSTS.ADMIN_EMAIL,
           receiver='orenavram@gmail.com', subject=f'MICROBIALIZER - A new job has been submitted: {run_number}',
           content=os.path.join(CONSTS.MICROBIALIZER_URL, 'results', run_number))

try:
    write_info_paragraph_to_html(output_html_path)

    content = ''
    content += f'{"#"*50}\n{ctime()}: A new CGI request has been recieved!\n'
    content += f'These are the keys that the CGI received:\n{"; ".join(sorted(form.keys()))}\n\n'
    content += 'Form values are:\n'
    for key in sorted(form.keys()):
        if 'data' not in key:
            content += f'{key} = {form[key]}\n'
        else:
            content += f'{key}'
    content += '\n'
    write_to_debug_file(cgi_debug_path, content)

    # extract form's values:
    user_email = form['email'].value.strip()
    example_page = form['example_page'].value

    job_title = ''
    if form['job_title'].value != '':
        job_title = form['job_title'].value.strip()

    write_to_debug_file(cgi_debug_path, f'Creating path {wd}')
    create_dir(wd)
    data_dir = os.path.join(wd, 'data')
    create_dir(data_dir)
    if example_page == 'no':
        # handling uploaded gzip file:
        original_file_name = form['data'].filename
        local_file_path = os.path.join(data_dir, original_file_name)
        with open(local_file_path, 'wb') as f:
            f.write(form['data'].value)
    else:
        write_to_debug_file(cgi_debug_path, f'{ctime()}: Copying example data...\n')
        # copy example data
        cp_cmd = f'cp {CONSTS.EXAMPLE_DATA} {data_dir}'
        write_to_debug_file(cgi_debug_path, f'Fetching: {cp_cmd}\n')
        os.system(cp_cmd)

    write_to_debug_file(cgi_debug_path, f'{ctime()}: ls of {wd} yields:\n{os.listdir(wd)}\n')

    # write_to_debug_file(cgi_debug_path, f'{ctime()}: Writing running parameters to html...')

    # write_running_parameters_to_html(output_html_path, job_title, actual_number_of_duplicates, urls_to_reads_files, files_names, MMU, len_threshold, qlty_threshold, number_of_clones_to_analyze, raw_data_suffix, mass_spec_seq, lib_file_name, f_umi, r_umi)

    # write_to_debug_file(cgi_debug_path, f'{ctime()}: Running parameters were written to html successfully.\n')

    # This is hidden field that only spammer bots might fill in...
    confirm_email_add = form['confirm_email'].value  # if it is contain a value it is a spammer.

    # parameters_file = os.path.join(wd, 'parameters.txt')
    # prepare_parameters_file(parameters_file, wd)

    output_dir = os.path.join(wd, 'output')
    create_dir(output_dir)

    cmds_file_path = os.path.join(wd, 'qsub.cmds')
    write_cmds_file(cmds_file_path, data_dir, output_dir, user_email, run_number)

    log_file_path = cmds_file_path.replace('cmds', 'log')
    # complex command with more than one operation (module load + python q_submitter.py)
    # submission_cmd = 'ssh bioseq@lecs2login "module load python/anaconda_python-3.6.4; python /bioseq/bioSequence_scripts_and_constants/q_submitter.py {} {} -q {} --verbose > {}"'.format(cmds_file_path, wd, queue_name, log_file_path)

    # simple command when using shebang header
    submission_cmd = f'ssh bioseq@lecs2login /bioseq/bioSequence_scripts_and_constants/q_submitter.py {cmds_file_path} {wd} -q bioseq --verbose > {log_file_path}'

    write_to_debug_file(cgi_debug_path, f'\nSSHing and SUBMITting the JOB to the QUEUE:\n{submission_cmd}\n')

    subprocess.call(submission_cmd, shell=True)

    if user_email != '':
        with open(os.path.join(wd, 'user_email.txt'), 'w') as f:
            f.write(f'{user_email}\n')

    write_to_debug_file(cgi_debug_path, f'\n\n{"#"*50}\nCGI finished running!\n{"#"*50}\n')

except Exception as e:
    msg = 'CGI crashed before the job was submitted :('
    with open(output_html_path) as f:
        html_content = f.read()
    html_content = html_content.replace('RUNNING', 'FAILED')
    html_content += '<br><br><br><center><h2><font color="red">' + msg + '</font><br><br>Please try to re-run your job or <a href="mailto:bioSequence@tauex.tau.ac.il?subject=MICROBIALIZER%20Run%20Number%2015249296875723">contact us</a> for further information</h2></center><br><br>\n</body>\n</html>\n'
    with open(output_html_path, 'w') as f:
        html_content = f.write(html_content)

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    write_to_debug_file(cgi_debug_path, f'\n{"$"*50}\n\n{msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\n{"$"*60}')

    # logger.info(f'Waiting {2*CONSTS.RELOAD_INTERVAL} seconds to remove html refreshing headers...')
    # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
    from time import sleep

    sleep(2 * CONSTS.RELOAD_INTERVAL)
    with open(output_html_path) as f:
        html_content = f.read()
    html_content = html_content.replace(CONSTS.RELOAD_TAGS, '')
    with open(output_html_path, 'w') as f:
        html_content = f.write(html_content)
