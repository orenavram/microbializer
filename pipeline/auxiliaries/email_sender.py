import argparse
from email.mime.text import MIMEText
from smtplib import SMTP


# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'


def send_email(logger, subject, content, email_addresses):
    if not email_addresses:
        logger.info('mail is empty, not sending')
        return

    if type(email_addresses) == str:
        email_addresses = [email_addresses]

    for email_address in email_addresses:
        try:
            send_email_internal(SMTP_SERVER, ADMIN_EMAIL, email_address, subject=subject, content=content)
            logger.info(f'sent email to {email_address} with subject {subject}')
        except:
            logger.exception(f'failed to send email to {email_address}')


def send_email_internal(smtp_server, sender, receiver, subject='', content=''):
    msg = MIMEText(content)
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = receiver
    s = SMTP(smtp_server)
    s.send_message(msg)
    s.quit()


if __name__ == '__main__':
    # This block will be executed only when you run it as your main program.
    # If this module is being imported from another script, this block won't be executed, however the function will be available...
    parser = argparse.ArgumentParser()
    parser.add_argument('smtp_server', help='SMTP server')
    parser.add_argument('sender', help='From whom the email will be sent')
    parser.add_argument('receiver', help='To whom the email will be sent')
    parser.add_argument('--subject', help='The subject of the email')
    parser.add_argument('--content', help='The content of the email', default='')
    args = parser.parse_args()
    send_email_internal(args.smtp_server, args.sender, args.receiver, args.subject, args.content)
