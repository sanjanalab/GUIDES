FROM ubuntu
MAINTAINER Joshua Meier <jmeier@college.harvard.edu>

RUN apt-get update
RUN apt-get install -y curl
RUN curl -sL https://deb.nodesource.com/setup | bash -
RUN apt-get install -y git nodejs libffi-dev libxml2-dev libxslt1-dev
RUN apt-get install -y python python-dev python-setuptools python-pip

RUN npm install -g less coffee-script

ADD . app

RUN cd app && pip install -r requirements.txt
RUN cd app && npm install

WORKDIR app
EXPOSE 5000
CMD python app.py
