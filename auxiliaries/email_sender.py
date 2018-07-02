def send_email(smtp_server, sender, receiver, subject='', content=''):
    from email.mime.text import MIMEText
    from smtplib import SMTP
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
    import logging, argparse
    logger = logging.getLogger('main') # use logger instead of printing
    parser = argparse.ArgumentParser()
    parser.add_argument('smtp_server', help='SMTP server')
    parser.add_argument('sender', help='From whom the email will be sent')
    parser.add_argument('receiver', help='To whom the email will be sent')
    parser.add_argument('--subject', help='The subject of the email')
    parser.add_argument('--content', help='The content of the email', default='')
    args = parser.parse_args()
    send_email(args.smtp_server, args.sender, args.receiver, args.subject, args.content)