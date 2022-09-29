#!/groups/pupko/modules/python-anaconda3.6.5/bin/python

import os
import shutil
import cgi
import cgitb
import subprocess
from time import time, ctime
from random import randint
import sys
sys.path.insert(0, '/bioseq/microbializer/pipeline/')  # for CONSTANTS and auxiliaries

import CONSTANTS as CONSTS
from auxiliaries.pipeline_auxiliaries import fix_illegal_chars_in_file_name
from auxiliaries.email_sender import send_email


def print_hello_world(output_path='', run_number='HAPPY FLOW'):

    hello_world_html = """
<html>
    <body>
        <h2>Hello World! """ + run_number + """</h2>
    </body>
</html>
    """
    if not output_path:
        print(hello_world_html)
    else:
        with open(output_path, 'w') as f:
            f.write(hello_world_html)


def write_to_debug_file(cgi_debug_path_f, content):
    cgi_debug_path_f.write(f'{ctime()}: {content}\n')


def write_html_prefix(output_path, run_number):
    with open(output_path, 'w') as output_path_f:
        output_path_f.write(f'''<html><head>
    
        <meta http-equiv="cache-control" content="no-cache, must-revalidate, post-check=0, pre-check=0" />
        <meta http-equiv="cache-control" content="max-age=0" />
        <meta http-equiv="expires" content="0" />
        <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
        <meta http-equiv="pragma" content="no-cache" />
        {CONSTS.RELOAD_TAGS}
    
        <title>{CONSTS.WEBSERVER_NAME} Job #{run_number}</title>
        <link rel="shortcut icon" type="image/x-icon" href="{CONSTS.WEBSERVER_URL}/pics/logo.gif" />
    
        <meta charset="utf-8">
        <!--<meta name="viewport" content="width=device-width, initial-scale=1">-->
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
        <link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">
    
        <link rel="stylesheet" href="{CONSTS.WEBSERVER_URL}/css/general.css">
        <link rel="stylesheet" href="../webpage/css/general.css">
    
        </head><body>
        <nav role="navigation" class="navbar navbar-inverse navbar-fixed-top">
            <div class="jumbotron" id="jumbo">
                <div class="container">            
                    <div class="row" id="title-row" align="center">
                        <div class="col-md-1">
                        </div>
                        <div class="col-md-10">
                            <span id="server-title">{CONSTS.WEBSERVER_NAME}</span>
                            <img src="{CONSTS.WEBSERVER_URL}/pics/logo.gif" id="nav_bar_image" class="img-rounded">
                            <br><span id="sub-title">{CONSTS.WEBSERVER_TITLE}</span>
                        </div>
                    </div>
                </div>       
            </div>
        </nav>
        <div id="behind-nav-bar-results">
        </div>
        <br><br><div class="container" style="font-size: 17px; {CONSTS.CONTAINER_STYLE}"  align="justify"><br> 
        <H1 align=center>Job Status: <FONT color='red'>\nQUEUED\n</FONT></h1>
        <br>
        {CONSTS.PROGRESS_BAR_TAG}
        {CONSTS.WEBSERVER_NAME} is now processing your request. This page will be automatically updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. A link to this page was sent to your email in case you wish to view these results at a later time without recalculating them. Please note that the results will be kept in the server for three months.
        <br><br></div>''')
        output_path_f.flush()


def write_running_parameters_to_html(output_path, identity_cutoff, e_value_cutoff, core_minimal_percentage, bootstrap, outgroup, job_title, original_file_name):
    with open(output_path, 'a') as output_path_f:
        output_path_f.write(f'<div class="container" style="{CONSTS.CONTAINER_STYLE}">')

        # optional param row
        if job_title != '':
            output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">')
            output_path_f.write(f'<b>Job title: </b>{job_title}')
            output_path_f.write('</div></div>')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">')
        output_path_f.write(f'<b>Archived data file name: </b>{original_file_name}')
        output_path_f.write('</div></div>')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">')
        output_path_f.write(f'<b>Maximal e-value cutoff: </b>{e_value_cutoff}')
        output_path_f.write('</div></div>')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">')
        output_path_f.write(f'<b>Identity minimal percent cutoff: </b>{identity_cutoff}%')
        output_path_f.write('</div></div>')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">')
        output_path_f.write(f'<b>Minimal percentage for core: </b>{core_minimal_percentage}%')
        output_path_f.write('</div></div>')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">')
        output_path_f.write(f'<b>Apply bootstrap: </b>{bootstrap.upper()}')
        output_path_f.write('</div></div>')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">')
        output_path_f.write(f'<b>Outgroup: </b>{outgroup if outgroup else "NONE"}')
        output_path_f.write('</div></div>')

        output_path_f.write('</div><br>')
        output_path_f.write('\n\n<!--result-->\n')
        output_path_f.flush()


def write_cmds_file(cmds_file, parameters, run_number):
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
    new_line_delimiter = '!@#'
    # the code contains features that are exclusive to Python3.6 (or higher)!
    required_modules = ' '.join(
        ['python/python-anaconda3.6.5-orenavr2'])
    with open(cmds_file, 'w') as f:
        f.write(f'module load {required_modules};')
        f.write(new_line_delimiter)
        f.write(f'{" ".join(["python", CONSTS.MAIN_SCRIPT, parameters])}\tmicrobializer_{run_number}')


def cleanup_is_running(cgi_debug_path_f, queues = ('pupkolab', 'pupkoweb')):
    for q in queues:
        write_to_debug_file(cgi_debug_path_f, f'\nfor q in queues\n')
        qstat = subprocess.Popen(('qstat', q), stdout=subprocess.PIPE)
        write_to_debug_file(cgi_debug_path_f, f'\nqstat\n')
        try:
            subprocess.check_output(('grep', 'cleanup_microb'), stdin=qstat.stdout).decode()
            write_to_debug_file(cgi_debug_path_f, f'\ncheck_output\n')
            # grep didn't crash, i.e., a cleanup is currently running.
            return True
        except subprocess.CalledProcessError:
            # grep returned nothing, i.e., no cleanup is currently running on this q.
            pass  # Keep checking other queues...
        except Exception as e:
            write_to_debug_file(cgi_debug_path_f, f'except Exception as e')
            write_to_debug_file(cgi_debug_path_f, f'\n{e.args[0]}\n')

    # none of the queues is running cleanup (none of them returned True)
    return False




def run_cgi():

    # prints detailed error report on BROWSER when backend crashes
    # This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
    cgitb.enable()

    # print_hello_world() # for debugging
    form = cgi.FieldStorage()  # extract POSTed object

    # adding 20 random digits to prevent users see data that are not their's
    run_number = str(round(time())) + str(randint(10 ** 19, 10 ** 20 - 1))
    # run_number = 'test'
    output_url = f'{CONSTS.WEBSERVER_RESULTS_URL}/{run_number}/{CONSTS.RESULT_WEBPAGE_NAME}'

    if form['example_page'].value == 'yes':
        run_number = 'example_data'
        output_url = os.path.join(f'{CONSTS.WEBSERVER_URL}/{run_number}/{CONSTS.RESULT_WEBPAGE_NAME}')

    # creating working directory
    wd = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, run_number)
    os.makedirs(wd, exist_ok=True)

    output_path = os.path.join(wd, CONSTS.RESULT_WEBPAGE_NAME)
    cgi_debug_path = os.path.join(wd, CONSTS.CGI_DEBUG_FILE_NAME)

    page_is_ready = os.path.exists(output_path)
    if not page_is_ready:
        write_html_prefix(output_path, run_number)  # html's prefix must be written BEFORE redirecting...

    print(f'Location: {output_url}')  # Redirects to the results url. MUST appear before any other print.
    print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
    sys.stdout.flush()  # must be flushed immediately!!!

    # email field should ALWAYS exist in the form (even if it's empty!!)
    # if it's not there, someone sent a request not via the website so they should be blocked.
    # confirm_email is hidden field that only spammer bots might fill in...
    if 'email' not in form or ('confirm_email' in form and form['confirm_email'].value != ''):
        exit()

    # Send me a notification email every time there's a new request
    send_email(smtp_server=CONSTS.SMTP_SERVER,
               sender=CONSTS.ADMIN_EMAIL,
               receiver=f'{CONSTS.OWNER_EMAIL}',
               subject=f'{CONSTS.WEBSERVER_NAME} - A new job has been submitted: {run_number}',
               content=f"{os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number, CONSTS.CGI_DEBUG_FILE_NAME)}\n"
                       f"{os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number, CONSTS.RESULT_WEBPAGE_NAME)}")

    try:
        cgi_debug_path_f = open(cgi_debug_path,'w')

        if 'email' in form and form['email'].value != '':
            write_to_debug_file(cgi_debug_path_f, f"{form['email'].value.strip()}\n\n")

        # for debugging
        write_to_debug_file(cgi_debug_path_f, f'{"#"*100}\n{ctime()}: A new CGI request has been received!\n')
        sorted_form_keys = sorted(form.keys())
        write_to_debug_file(cgi_debug_path_f, f'These are the keys that the CGI received:\n{"; ".join(sorted_form_keys)}\n\n')
        write_to_debug_file(cgi_debug_path_f, 'Form values are:\n')
        for key in sorted_form_keys:
            if 'data' not in key:
                write_to_debug_file(cgi_debug_path_f, f'{key} = {form[key]}\n')
        for key in sorted_form_keys:
            if 'data' in key:
                write_to_debug_file(cgi_debug_path_f, f'100 first characters of {key} = {form[key].value[:100]}\n')

        # extract form's values:
        email = form['email'].value.strip()
        job_title = form['job_title'].value.strip()
        identity_cutoff = form['identity_cutoff'].value.strip()
        e_value_cutoff = form['e_value_cutoff'].value.strip()
        core_minimal_percentage = form['core_minimal_percentage'].value.strip()
        bootstrap = form['bootstrap'].value.strip()
        outgroup = form['outgroup'].value.strip()

        if form['example_page'].value == 'no':
            write_to_debug_file(cgi_debug_path_f, f'\n{"#"*80}\nuploading data\n')

            file_name = form['data'].filename
            write_to_debug_file(cgi_debug_path_f, f'file name is:\n{file_name}')

            if not file_name.endswith(('zip', 'tar.gz', 'gz')):
                write_to_debug_file(cgi_debug_path_f, f'FILE FORMAT IS ILLEGAL!!')
                raise ValueError

            original_file_name = file_name
            file_name = fix_illegal_chars_in_file_name(file_name)

            data = form['data'].value
            write_to_debug_file(cgi_debug_path_f, f'{file_name} first 100 chars are: {data[:100]}\n')

            data_path = os.path.join(wd, file_name)
            with open(data_path, 'wb') as data_f:
                data_f.write(data)
            write_to_debug_file(cgi_debug_path_f, f'Uploaded data was saved to disk successfully\n')
        else:  # example data
            original_file_name = file_name = CONSTS.EXAMPLE_DATA_FILE_NAME
            data_path = os.path.join(wd, file_name)
            write_to_debug_file(cgi_debug_path_f, f'Copying example data FROM {CONSTS.EXAMPLE_DATA} TO {data_path}\n')
            try:
                shutil.copy(CONSTS.EXAMPLE_DATA, data_path)
                write_to_debug_file(cgi_debug_path_f, f'File was copied successfully to {data_path}\n\n')
            except IOError as e:
                write_to_debug_file(cgi_debug_path_f, f'{"#" * 50}\nFailed to copy for example data due to the following reason:\n{e.args[0]}\n{"#" * 50}\n')
                raise e

        write_to_debug_file(cgi_debug_path_f, f'ls of {wd} yields:\n{os.listdir(wd)}\n')

        write_to_debug_file(cgi_debug_path_f, f'{ctime()}: write_running_parameters_to_html...\n')

        if not page_is_ready:
            write_running_parameters_to_html(output_path, identity_cutoff, e_value_cutoff, core_minimal_percentage, bootstrap, outgroup, job_title, original_file_name)

        write_to_debug_file(cgi_debug_path_f, f'{ctime()}: Running parameters were written to html successfully.\n')

        # send main jobs ONLY TO ONE MACHINE so there won't be a dead lock
        # (cluster might be full with main jobs waiting for sub jobs but they are in qw mode...)
        queue_name = 'pupkor'  # all pupko machines on power
        queue_name_for_subjobs = 'pupkowebr'
        # queue_name_for_subjobs = '"pupkotmpr -l nodes=compute-0-18"'
        # queue_name = '"pupkowebr -l nodes=compute-0-292"'  # run main thread on pupkolab machine to prevent deadlock
        # queue_name_for_subjobs = '"pupkowebr -p -1"'  # reduced priority

        parameters = f'"{data_path}" ' \
                     f'{os.path.join(wd, "outputs")} ' \
                     f'--identity_cutoff {identity_cutoff} ' \
                     f'--e_value_cutoff {e_value_cutoff} ' \
                     f'--core_minimal_percentage {core_minimal_percentage} ' \
                     f'--bootstrap {bootstrap} ' \
                     f'--outgroup "{outgroup}" ' \
                     f'--queue_name {queue_name_for_subjobs} '
        if email != '':
            parameters += ' ' + f'--email {email}'

        cmds_file = os.path.join(wd, 'qsub.cmds')
        write_cmds_file(cmds_file, parameters, run_number)

        job_id_file = os.path.join(wd, 'job_id.txt')
        # simple command when using shebang header
        submission_cmd = f'{CONSTS.Q_SUBMITTER_PATH} {cmds_file} {wd} -q {queue_name} --verbose > {job_id_file}'

        write_to_debug_file(cgi_debug_path_f, f'\nSUBMITTING JOB TO QUEUE:\n{submission_cmd}\n')
        if not page_is_ready:
            subprocess.call(submission_cmd, shell=True)

        if email != '':
            with open(os.path.join(wd, CONSTS.EMAIL_FILE_NAME), 'a') as email_f:
                email_f.write(f'{email}\n')

            job_name = f"Job title: {job_title}\n" if job_title else ''
            notification_content = f"Your submission details are:\n\n{job_name}" \
                                   f"Dataset name: {file_name}\n" \
                                   f"Maximal e-value cutoff: {e_value_cutoff}\n" \
                                   f"Identity minimal percent cutoff: {identity_cutoff}%\n" \
                                   f"Minimal percentage for core: {core_minimal_percentage}%\n" \
                                   f"Apply bootstrap procedure: {bootstrap.upper()}\n" \
                                   f"Outgroup: {outgroup}\n\n" \
                                   f"Once the analysis will be ready, we will let you know! " \
                                   f"Meanwhile, you can track the progress of your job at:\n{os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number, CONSTS.RESULT_WEBPAGE_NAME)}\n\n"

            # Send the user a notification email for their submission
            try:
                send_email(smtp_server=CONSTS.SMTP_SERVER,
                           sender=CONSTS.ADMIN_EMAIL,
                           receiver=f'{email}',
                           subject=f'{CONSTS.WEBSERVER_NAME} - Your job has been submitted! (Run number: {run_number})',
                           content=notification_content)
            except:
                write_to_debug_file(cgi_debug_path_f, f'\nFailed sending notification to {email}\n')

        write_to_debug_file(cgi_debug_path_f, f'\n\n{"#"*50}\nCGI finished running!\n{"#"*50}\n')

        # run periodical cleanup on the background...
        try:
            # if not cleanup_is_running(cgi_debug_path_f):
            submission_cmd = f'{CONSTS.Q_SUBMITTER_PATH} {CONSTS.WEBSERVER_RESULTS_DIR}/cleanup_microbializer.cmds {CONSTS.WEBSERVER_RESULTS_DIR}'
            write_to_debug_file(cgi_debug_path_f, f'\nSUBMITTING CLEANUP:\n{submission_cmd}\n')
            subprocess.call(submission_cmd, shell=True)
        except Exception as e:
            write_to_debug_file(cgi_debug_path_f, f'\n\n{"%"*50}\n'
                                                  f'Crashed when tried to fetch folders cleanup for some reason...\n'
                                                  f'Error message: {e.args}')

        cgi_debug_path_f.close()

    except Exception as e:
        msg = f'{CONSTS.WEBSERVER_NAME} could not save your data 😞<br>\n' \
              'Most likely because your data contain non <a href="http://www.asciitable.com/" target="_blank">ASCII characters</a>.'
        with open(output_path) as f:
            html_content = f.read()
        html_content = html_content.replace('QUEUED', 'FAILED')
        html_content = html_content.replace(' active', '')
        html_content += f'<br><br><br>' \
                        f'<div class="container" style="{CONSTS.CONTAINER_STYLE}"><h3>' \
                        f'<font color="red">{msg}</font><br><br>' \
                        f'Please try to re-run your job in a few minutes or <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME}%20Run%20Number%20{run_number}">contact us</a> for further information' \
                        f'</h3></div><br><br>\n</body>\n</html>\n'
        with open(output_path, 'w') as f:
            f.write(html_content)

        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        with open(cgi_debug_path, 'w') as cgi_debug_path_f:
            write_to_debug_file(cgi_debug_path_f, f'\n{"$"*100}\n\n{msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\ne.args[0]: {e.args[0]}\n\n{"$"*100}')

        # Send me a notification email every time there's a failure
        try:
            email = form['email'].value.strip() if form['email'].value.strip() else 'NO_EMAIL'
        except:
            email = 'NO_EMAIL'
        send_email(smtp_server=CONSTS.SMTP_SERVER,
                   sender=CONSTS.ADMIN_EMAIL,
                   receiver=f'{CONSTS.OWNER_EMAIL}',
                   subject=f'{CONSTS.WEBSERVER_NAME} job {run_number} by {email} has been failed!',
                   content=f"{email}\n\n{os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number, CONSTS.RESULT_WEBPAGE_NAME)}\n"
                           f"\n{os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number, CONSTS.CGI_DEBUG_FILE_NAME)}")

        # logger.info(f'Waiting {2*CONSTS.RELOAD_INTERVAL} seconds to remove html refreshing headers...')
        # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
        from time import sleep

        sleep(2 * CONSTS.RELOAD_INTERVAL)
        with open(output_path) as f:
            html_content = f.read()
        html_content = html_content.replace(CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')
        with open(output_path, 'w') as f:
            f.write(html_content)

    # logging submission
    with open(CONSTS.SUBMISSIONS_LOG, 'a') as f:
        f.write(f'{email}\t{run_number}\t{ctime()}\n')

    with open(cgi_debug_path, 'a') as f:  # for cgi debugging
        f.write(f'{ctime()}: Submission was documented in \n')


if __name__ == '__main__':
    run_cgi()
