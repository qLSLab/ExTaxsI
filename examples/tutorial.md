ExTaxsI Tutorial
=======

This is the tutorial for ExTaxsI tool. The procedures described were run with the files in the examples directory of our Github page. \
In particular, the tutorial will explain the steps of installation, manual query investigation, query through a file (list of accessions, list of taxids
or list of organisms name) and plotting.
Please, for any doubt or issue, contact us through github issue messages or sending an e-mail to g.agostinetto@campus.unimib.it.

Before starting
---------------------
1) Download Python 3 if you haven’t already at https://www.python.org/
2) Clone github repository (link: https://github.com/qLSLab/ExTaxsI) and:
* Download ZIP from ExTaxsI github home page, then decompress \
or 
* Use git from command line:
git clone https://github.com/qLSLab/ExTaxsI

3) Create the conda environment with the script below, modifying 'myExTaxsIenv' with the name that you prefer: \
conda create --name myExTaxsIenv --file requirements.txt --channel default --channel etetoolkit --channel plotly

4) Modify settings.ini at:
* Line 8: insert your entrez e-mail
* Line 17: insert your entrez API key

To create entrez email and API key, you must create a NCBI account (follow the tutorial here https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

Warning settings
---------------------
* If you have issues with HTTP connection, make sure to fill the NCBI keys in settings.ini
* This tutorial is built only to demonstrate ExTaxsI functioning: results are not discussed and can change if you consider different versions of the Taxonomy DB or different moments for the download.

Run ExTaxsI
---------------------
Run the script below to run the tool:
$ python ExTaxsI.py 

Suggestions: 
ExTaxsI implemented the tab compiler, you can use it to insert input files from other directories (e.g. >>> download/A_accession_list_example.tsv)
do not use spaces ‘ ‘ when you give names to file (we suggest to use underscores ‘_’).

The tool will ask if you want to update the Taxonomy DB, press Y or N to answer.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_0a.png)

Then the tool will plot the main menu with the three main modules: Database, Visualization and Taxonomy ID Converter.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_0b.png)

You can press the module that you prefer. \
**To follow the tutorial, press 1.**

Database module
---------------------
This module allows to create from a NCBI query or from input files the following outputs:
* sequences in fasta format
* accessions + taxonomy in tsv format
* dataframe with genes found
* a PNG image with the top 10 most abundant genes

Pressing 1 from main menu, you will see the following question:
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1a.png)

Press m for manual query and press enter. \
You will see the following output:
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1b.png)

