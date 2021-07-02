ExTaxsI step-by-step tutorial and example usage
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

3) We suggest to install Conda or Mininconda to run ExtaxsI - here you can find how to install Miniconda3: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

4) Once Miniconda is installed, create the conda environment with the script below, modifying 'myExTaxsIenv' with the name that you prefer: \
conda create --name myExTaxsIenv --file requirements.txt --channel default --channel etetoolkit --channel plotly

5) Modify settings.ini at:
* Line 8: insert your entrez e-mail
* Line 17: insert your entrez API key

To create entrez email and API key, you must create a NCBI account (follow the tutorial here https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

Warning settings
---------------------
* If you have issues with HTTP connection, make sure to fill the NCBI keys in settings.ini
* This tutorial is built only to demonstrate ExTaxsI functioning: results are not discussed and can change if you consider different versions of the Taxonomy DB or different moments for the download.

Input files explanation
---------------------
* A_accession_list_example.tsv: list of accessions (one accession per line), the file should be in a csv or tsv format, all the accessions will be in the first column;
* T_taxid_list_example.tsv: list of TaxIDs (one TaxID per line), the file should be in a csv or tsv format, all the accessions will be in the first column;
* O_organism_list_example.tsv: list of organism names (one organism name per line), the file should be in a csv or tsv format, all the accessions will be in the first column.

NCBI query explanation
---------------------
Nucleotide + Entrez query = https://www.ncbi.nlm.nih.gov/homologene/query-tips/

Run ExTaxsI
---------------------
Activate the conda environment with:
$ conda activate myExTaxsIenv

Run the script below to run the tool:
$ python ExTaxsI.py 

Suggestions: 
ExTaxsI implemented the tab compiler, you can use it to insert input files from other directories (e.g. >>> download/A_accession_list_example.tsv)
do not use spaces ‘ ‘ when you give names to file (we suggest to use underscores ‘_’).

The tool will ask if you want to update the Taxonomy DB, press Y or N to answer.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_0a.png)

Then the tool will plot the main menu with the three main modules: Database, Visualization and Taxonomy ID Converter.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_0b.png)

You can press the module that you prefer.
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

### Manual query
Press m for manual query and press enter. \
You will see the following output:
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1b.png)

You can see example queries. Type ‘carabus convexus’ to follow this tutorial.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1c.png)

\

As you can see, Database module has several options.
Enter the number for the option that you prefer, all the outputs will be saved in 'download' directory.
Press 4 to obtain all the possible outputs.

Pressing **the option 3** ExTaxsI will create a dataframe with all the genes available for the input query. If you want, you can create a PNG image with the top 10 most abundant genes.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1d.png)

### Query through a file
You can submit to ExTaxsI your personal lists as accessions, TaxIDs or organism name. In the example directory you can find example files (with the correct format).
In particular, you have to create a list of accessions, TaxIDs or organism name and save as csv or tsv format.
The example files are the following:
* A_accession_list_example.tsv
* T_taxid_list_example.tsv
* O_organism_list_example.tsv

You can press A for a list of accession, T for a list of TaxIDs and O for a list of organism name.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1e.png)

Rembember to type the directory in which you saved the input file.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1f.png)

After giving the input file name with the directory in which it is saved, ExTaxsI will ask which type of output you want to download 
(please see the previous images for a description of the options available).
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1g.png)

Choosing **T option**, ExTaxsI will ask if you want a set of genes - press 0 if you want all.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1h.png)

As for A option, you must enter the file (with the directory) and choose the options form the menu.
The same will happen with the O option.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_1i.png)

Visualization module
---------------------
The Visualization module (option 2 from Main Menu) allows to create plots to visualize the taxonomy distribution of your query or your input files. In particular, it can generate:
* scatter plot
* world map plot
* sunburst plot
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_2a.png)

### Scatter plot
You can generate a scatterplot with the following options: a taxonomy file created with ExTaxsI (from download directory), doing a query or with an external tsv or csv file (see examples for details).
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_2b.png)

Using option 1 and 3, you can use as input a tsv or csv file with accessions and taxonomy (see the example directory), as we described into the **Database module**.
**INSERIRE IN EXAMPLE IL FILE** examples/ace2_taxid7742_taxonomy.tsv
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_2c.png)
With option 2, you can digit your query (like in Database module, see the previous slides for details). Then ExTaxsI will create the file with accessions and taxonomy and it will plot it.

For each option, it is possible to **filter the number of accessions for each species** to be considered to generate the plot. Extaxsi will ask using the _filter value_.

![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/txid7742_ace2_taxonomy_scatterplot.png)

### World map plot
You can generate a world map plot with two modalities: manual query or through a file, see Database module tutorial for details (example files for input were provided in the tutorial folder).
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_2e.png)
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_2f.png)

INSERIRE WORLD MAP PLOT DI ESEMPIO

### Sun burst plot
The modalities to create the sunburst plot are the same of scatterplot.
Please for details refers to scatterplot examples tutorial.

INSERIRE SUNBURST PLOT DI ESEMPIO

Taxonomy ID Converted module
---------------------
The module 3 allows to convert your TaxIDs into a 6-rank taxonomy or obtained from an organism name the corresponding TaxID. The query can be done both with manual and input file query (see example directory for input files used in the tutorial).
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_3a.png)

### Manual query
Choosing manual query, you can then decide if you want to convert TaxID to Taxonomy name (6-rank) or viceversa.
**Warning**: If you submit the TaxID of higher rank than species (e.g. order) the following rank will be filled as NA (Not Assigned).

![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_3b.png)

If you want to submit multiple TaxIDs use ‘+’ as a separator.
Then ExTaxsI will return the complete taxonomy.

![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_3c.png)

Choose option 2 to convert organism name to TaxID.
If you want to submit multiple organisms use ‘+’ as a separator (space only between genus and species).
Then ExTaxsI will return the corresponding TaxIDs.

![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_3d.png)

### Query through an input file
Press f for input file query. In particular, you can use a list of TaxIDs or a list of organism names (please see example directory for example input files).
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_3e.png)
Press T for a list of TaxIDs or O for a list of organism names (file must have tsv or csv extension).

The tool asks if there is an header into the file, then it saves the results into the download directory.
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/tutorial_3f.png)
