# Graphical User Interface for DNA Editing Screens (GUIDES)

**Description**:  Genome-scale CRISPR-Cas9 knockout libraries have emerged as popular tools for unbiased, phenotypic screening but easy-to-use software for designing customized guide RNA libraries has lagged behind. GUIDES is a web application for automated, tissue-specific library design. It integrates on-target scores, expression data, and protein structure into an easily accessible web interface.

  - **Technology stack**: Python/Flask backend with Angular.js front-end.
  - **Status**:  Beta (1.0)
  - **Links to production or demo instances**: Hosted on Zhang Lab's [Genome Engineering](guides.genome-engineering.org) website.
  - Recently, library design tools such as CRISPR Library Designer have used on-target efficiency scoring to help automate guide design. In addition to on-target scores, our tool integrates expression data and protein structure, while also providing a more easily accessible web interface.


**Screenshots**:

![](https://raw.githubusercontent.com/fengzhanglab/GUIDES/master/Screenshot.png)


## Dependencies

Most of the front-end code is written in Coffeescript.

Redis/Celery are used for local storage.

Static data files which are too large for hosting on GitHub must also be added. Installation is described below.

## Installation

After pulling this repository, copy `data/` from the [Zhang Lab Dropbox]() and move it to `CRISPR_Library_Design/static`.

Alternatively, run [./install.sh](install.sh).

## Usage

Next, open 3 terminal windows. In each window, run the following commands:

```bash
source venv/bin/activate
export DEV=false
export PORT=YOUR_PORT
export SMTPserver=YOUR_SMTP_SERVER
export SMTPsender=donotreply@YOUR_DOMAIN
export SMTPpassword=YOUR_SMTP_PASSWORD
```

Then, run the following 3 commands, in this order (one in each window):

1. `./run-redis.sh`
2. `celery worker -A app.celery -f celery_log.txt --detach --autoreload`
3. `./run.sh`

Detailed instructions on how to install, configure, and get the project running.
This should be frequently tested to ensure reliability. Alternatively, link to
a separate [INSTALL](INSTALL.md) document.

## How to test the software

Code for running the tests in our [manuscript]() is located in `static/data/functional_tests`. Some data analysis is required to produce the plots included in our paper.

## Known issues

None currently.

## Getting help

Feel free to reach out to the Zhang Lab team at [cld@mit.edu](cld@mit.edu).

**Example**

Here is an [example library]() targeting 500 genes involved in chromatin regulation with 6 guides per gene.

## Getting involved

Please [get in touch](cld@mit.edu) if you are interested in contributing, so we can coordinate our efforts. In particular, we are interested in expanding our tool to include more species.

----

## Open source licensing info
Copyright (c) 2015-2016 Joshua Meier. Released under AGPLv3. See LICENSE.txt for details.

If you're interested in getting access to this system under a different license, please contact me.

----
