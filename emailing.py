# send emails
import requests
import os

def send_completed_run(address, link):

  key = os.environ.get('mailgun_key')
  sandbox = os.environ.get('mailgun_sandbox')
  if address:
    recipient = address + ", cld.history@gmail.com"
  else:
    recipient = "cld.history@gmail.com"

  content = """Hello {0},
  GUIDES (Graphical User Interface for DNA Editing Screens) has completed processing your library.
  Results are available here: http://guides.sanjanalab.org/#/designer/{1}

  Please contact us at guides@mit.edu if you have any questions or feedback.

  The GUIDES Team
  """.format(address, link)

  request_url = 'https://api.mailgun.net/v3/{0}/messages'.format(sandbox)
  request = requests.post(request_url, auth=('api', key), data={
      'from': 'donotreply@cld.genome-engineering.org',
      'to': recipient,
      'subject': "GUIDES Run Completed",
      'text': content
  })

# import sys
# import os
# import re
# from smtplib import SMTP
# from email.mime.text import MIMEText

# def send_completed_run(address, link):
#   SMTPserver = os.environ.get('SMTPserver')
#   sender = os.environ.get('SMTPsender')
#   password = os.environ.get('SMTPpassword')
#   destination = [address, "cld.history@gmail.com"]

#   text_subtype = 'plain'
#   subject = "GUIDES Run Completed"

#   content = """Hello {0},
#   GUIDES (Graphical User Interface for DNA Editing Screens) has completed processing your library.
#   Results are available here: http://guides.sanjanalab.org/#/designer/{1}

#   Please contact us at guides@mit.edu if you have any questions or feedback.

#   The GUIDES Team
#   """.format(address, link)

#   try:
#     msg = MIMEText(content, text_subtype)
#     msg['Subject'] = subject
#     msg['From'] = sender
#     msg.add_header('reply-to', 'guides@mit.edu')

#     conn = SMTP(SMTPserver)
#     conn.set_debuglevel(False)
#     conn.login(sender, password)
#     try:
#       conn.sendmail(sender, destination, msg.as_string())
#     finally:
#       conn.quit()

#   except Exception, exc:
#     print "mail failed; %s" % str(exc)
