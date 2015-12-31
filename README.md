# CRISPR Library Designer

A new workflow for the custom design of CRISPR Libraries.

To run the server locally, follow these steps.

1) Execute
git clone https://github.com/joshim5/CRISPR-Library-Designer.git

2a) Download and unzip the two tarballs from my Dropbox directory, https://www.dropbox.com/sh/pfn4f60wg13sjt5/AADb57WJdkMbWzNZR3Qu_C07a?dl=0. 

2b) Make sure the resulting directories are called gtex/ and pre_processed/.

2c) Move gtex/ and pre_processed/ to CRISPR-Library-Designer/static/data/.

3) Download the GRCh37 sequences from ensembl. Put the tarballs inside CRISPR-Library-Designer/static/data/GRCh37/.

4) In termainl, naviate to CRISPR-Library-Designer/. Execute
pip install -r requirements.txt

5) Execute these commands:
export PORT=8000
export DEV=true

6)Execute
./run.sh

Fainlly, open localhost:8000 in your browser and you'll see the library designer :)

- Joshua Meier, Zhang Lab, 2015