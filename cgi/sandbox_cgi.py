#!/data/shared/python/anaconda3-5.1.0/bin/python3.6

# old shebang: #!/shared/python/anaconda3.5/bin/python

#THIS IS A SANDBOX!!!!!!
#THIS IS A SANDBOX!!!!!!
#THIS IS A SANDBOX!!!!!!

# Common error: "[Sat Apr 28 11:22:27.659177 2018] [http:error] [pid 29692] [client 132.66.1.18:21082]
# AH02429: Response header name '<!--' contains invalid characters, aborting request, referer: http://localhost:63343/webpage/index.html"
# but it's not always the real reason!!
# Solution: It may happen when trying to get undefined form key!! For more details see:
# https://stackoverflow.com/questions/45375693/apache2-response-header-name-contains-invalid-characters-aborting-reque

# Common error: "[Sun Jun 03 13:18:14.345457 2018] [cgi:error] [pid 22240] [client 132.66.1.18:61229] End of script output before headers: MICROBIALIZER_cgi.py"
# Solution: make sure there are execution permissions!! chmod 744
# https://www.1and1.com/cloud-community/learn/web-server/server-management/how-to-fix-http-error-code-500-internal-server-error/

import cgi
import cgitb
import os
import sys

#sys.path.append('/bioseq/bioSequence_scripts_and_constants')
sys.path.append('/bioseq/MICROBIALIZER/auxiliaries')
import CONSTANTS as CONSTS

def hello_world(output_path='', run_number='NO_RUN_NUMBER'):
    with open(output_path, 'w') as f:
        f.write("""<html><body><h1 align=center><FONT color='red'>""" + run_number + """ was POSTed successfully!!!</FONT></h1>""")

def write_running_parameters_to_html(output_path, form):
    with open(output_path, 'a') as f:

        f.write("""<font face=Verdana><u><h4>Running Parameters:</h4></u></font>""")
        f.write(f'These are the keys that the CGI received:<br>{"; ".join(sorted(form.keys()))}<br><br>')
        for key in sorted(form.keys()):
            if 'run' not in key:
                f.write(f'{key} = {form[key]}<br>')
        for key in sorted(form.keys()):
            if 'run' in key:
                f.write(f'100 first characters of {key} = ')
                f.write(f'{form[key].value[:100]}<br>')
        f.write('<br><br>')
        f.write('</body></html>')

# def print_running_parameters_to_html(form):
#     print("""<font face=Verdana><u><h4>Running Parameters:</h4></u></font>""")
#     for key in form:
#         if 'run' not in key:
#             print('<ul>{} = {}<br></ul>'.format(key, form[key]))
#     print('<br><br>')
#     print('</body></html>')
#
# def print_hello_world(run_number='NO_RUN_NUMBER'):
#     print("""<html><body><h1 align=center><FONT color='red'>""" + run_number + """ was POSTed successfully!!!</FONT></h1>""")


# prints detailed error report on BROWSER when cgi crashes
# This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
cgitb.enable()

# print_hello_world() # for debugging
form = cgi.FieldStorage()  # extract POSTed object

run_number = 'sandbox'

results_url = os.path.join(CONSTS.MICROBIALIZER_RESULTS_URL, run_number)
output_url = os.path.join(results_url, 'output.html')

wd = os.path.join(CONSTS.MICROBIALIZER_RESULTS_DIR, run_number)
os.makedirs(wd, exist_ok=True)

output_path = os.path.join(wd, 'output.html')
cgi_debug_path = os.path.join(wd, 'cgi_debug.txt')

hello_world(output_path, run_number)
write_running_parameters_to_html(output_path, form) # html's prefix must be written BEFORE redirecting...

print('Location: ' + output_url)  # Redirects to the results url. MUST appear before any other print.
print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
sys.stdout.flush()  # must be flushed immediately!!!

# print_hello_world(output_html_path, run_number)
# print_running_parameters_to_html(form) # html's prefix must be written BEFORE redirecting...





