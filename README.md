ExTaxsI
=======

![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/Project%20Exta.png)

Project overview
----------------
ExTaxsI is a bioinformatic tool aimed to elaborate and visualize molecular and taxonomic informations.
This open-source user friendly project, written in Python 3.7, allows the creation of interactive plots starting from NCBI search query or directly from offline taxonomic files.

ExTaxsI has multiple functions:

* DATABASE: creation of multi FASTA files composed by nucleotide sequences, taxonomic lists, genes and accessions, starting from manual inputs or csv/tsv files.

* VISUALIZATION: creation of interactive plots, such as scatter plot, sun burst plot and world map, starting from DATABASE output or external sources.

* ID CONVERSION: conversion of TaxID into 6-ranks taxonomy and vice versa; it can convert single manual inputs, takes multiple inputs joined together by a plus sign or tsv/csv file with a list of taxIDs.

Hardware requirements
---------------------
Minimum hardware requirements:
no specific requirements are needed for ExTaxsI installation, however for the correct functioning of the software we suggest the following:

* RAM: 4GB
* CPU: quad-core or more.

Installation instructions in a nutshell
-------------------------
**1- Download Python 3:**

https://www.python.org/

**2- INSTALL EXTERNAL LIBRARIES:**

To install external libraries open a terminal (prompt for windows users) in the extaxsi folder or navigate to the folder with the following:
* cd path/to/ExTaxsI-folder
* cd path\to\ExTaxsI-folder (windows-users)

Now that you are inside the ExTaxsI-folder, if you are a conda user (best option), run the following command:
* conda create --name myExTaxsIenv --file requirements.txt --channel default --channel etetoolkit --channel plotly

Otherwise (less recommended), run the following command:
* pip3 install -r pip_requirements.txt

It is also required to install plotly-orca using the following instructions:
https://github.com/plotly/orca

**3- CUSTOMIZE YOUR SETTINGS:**

Before starting to use ExTaxsI, the settings.ini file should be customized:
* entrez_email: insert your email;
* api_key: insert your api_key created from NCBI;

In order to not overload the NCBI servers, by entering your API key in the setting file, NCBI admits maximum 10 requests/second for all activities from that key.

Here is the reference: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

How to run ExTaxsI
------------------
Open a terminal and go to the ExTaxsI folder directory by running the entire path:

linux:

* cd path/to/ExTaxsI-folder

windows (prompt):

* cd path\to\ExTaxsI-folder

mac:

* cd path/to/ExTaxsI-folder

Now that you are in the right directory, run the following command to start:

* python ExTaxsI.py

Please visit the tutorial directory of this Github page for a step-by-step installation and example usage.

Operating instructions
----------------------

The first time that you run ExTaxsI, the program will take time to download the local database which you can update as needed on startup.

**Which module do you want to use?**

Choose the module you’re interested in by entering the correlated number:

1. Database creation module: taxonomy and FASTA files download;
2. Visualization module: scatter plot and world map plot from taxonomy files or queries;
3. Taxonomy IDs converter: conversion of taxID to 6-ranks taxonomy and vice versa;


**Module 1: Database creation module;
taxonomy and fasta files download.**

When organism name list, Ids or accessions are less than 2500 the search key is composed by a single query, otherwise query will be splitted in groups of 2500 generating temporary files, which would be deleted at the end of the process.

Output file (standard format) will be saved in Download folder.

Available formats:
* multi-FASTA file format (NCBI standard format, with header followed by nucleotide sequence, accession, name, code, gene)
* TSV format


**Module 2: Statistical module; ScatterPlot and world map from taxonomy files or queries**

Module2’s required data, can be uploaded in several ways:
* manually, by entering a query;
* uploading module1’s taxonomy output file;
* uploading file from external sources containing taxonomy lists.

It's possible to do 3 types of interactive plots:
* scatter plot: uses taxonomy as input to produce a graph that indicates the quantity of each individual taxonomic unit and which taxa are present;
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/aves%20scatterplot%20COX1.png)
* sunburst plot: uses taxonomy but creates an expansion pie that allows to explore taxonomy in depth with less weight on the quantity of each individual taxonomic unit;
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/sunburst%20odonata.png)
* world map plot: uses the country metadata of accessions data to produce a map indicating the position of each taxon found;
![alt text](https://github.com/qLSLab/extaxsi/blob/master/images/worldmap.png)  

Output file format: html

Output file folder: download

**Module3:  Taxonomy ID converter.**

You can convert taxonomy ids from a file or by manual input into full 6 main ranks taxonomy (phylum;class;order;family;genus;species).

__Note__: list values must be placed within the first column.

Copyright and licensing information
-----------------------------------

Copyright (c)

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Contacts
-----------------------------------------------------
* Name: __Giulia Agostinetto__
* e-mail: __g.agostinetto@campus.unimib.it__
* Name: __Alberto Brusati__
* e-mail: __alberto.brusati@gmail.com__
* Name: __Anna Sandionigi__
* e-mail: __anna.sandionigi@unimib.it__

Credits and acknowledgments
---------------------------

Throughout the creation of this tool, we relied heavily on contributions from professors and college students. Their input was invaluable, and we want to take a moment to thank them and recognize them for all of their hard work:
* __Adam Chahed__, Ex-student at University of Milano-Bicocca.
* __Elena Parladori__, Ex-student at University of Milano-Bicocca.
* __Bachir Balech__, Researcher at CNR IBIOM.
* __Dario Pescini__, Associate Professor at University of Milano-Bicocca.
