EXTAXSI
=======

![alt text](https://github.com/albertobrusati/extaxsi_alberto/blob/master/Project%20Exta.png)

Project overview
----------------
ExTaxsI is a bioinformatic tool aimed to elaborate and visualize molecular and taxonomic informations.
This open-source user friendly project, written in Python 3.7, allows the creation of interactive plots starting from NCBI search query or directly from offline taxonomic files.

ExTaxsI has multiple functions:

* DATABASE: creation of multi FASTA files composed by nucleotide sequences, taxonomic lists, genes and accessions, starting from manual inputs or csv/tsv files.

* VISUALIZATION: creation of interactive plots, such as scatterplot, sunburst and worldmap, starting from dbs output or external sources.

* ID CONVERSION: conversion of taxID to 6ranks taxonomy and vice versa; can convert single manual inputs, take multiple inputs joined together by a plus sign or tsv/csv file with a list of taxIDs.

Hardware requirements
---------------------
Minimum hardware requirements:
no specific requirements are needed for ExTaxsI installation, otherwise for the correct functioning of the software we suggest the following:

* RAM: 4GB
* CPU: quad-core or more.

Installation instructions
-------------------------
**1- Download Python:**

https://www.python.org/

**2- INSTALL EXTERNAL LIBRARIES:**

To install external libraries open terminal;

if you are a conda user, run the following command:
* conda create - -name myenv - -file conda_req.txt

otherwise, run the following command:
* pip install -r requirements.txt

**3- CUSTOMIZE YOUR SETTINGS:**

Before starting to use ExTaxsI, the settings.ini file should be customized:
* entrez_email: insert your email;
* api_key: insert your api_key created from NCBI;

In order to not overload servers, by entering your APY key in the setting file, NCBI admits maximum 10 requests/second for all activity from that key;


How to run extaxsi
------------------
Open terminal to the ExTaxsI folder directory by running the entire path:

linux:

* cd /path/to/ ExTaxsI folder

windows:

* cd /d D: \ ExTaxsI folder

mac:

* cd ~/path/to/ ExTaxsI folder

Now that you are in the right directory, run the following command to start:

* python extaxsi.py

Operating instructions
----------------------

The first time that you run ExTaxsI, the program will take time to download the local database which you can update as needed on startup.

**Which module do you want to use?**

Choose the module you’re interested in by entering the correlated number:

1. Database creation module: taxonomy and FASTA files download;
2. Statistical module: scatter plot and world map from taxonomy files or queries;
3. Taxonomy IDs converter: conversion of taxID to 6ranks taxonomy and vice versa;


**Module 1: Database creation module;
taxonomy and fasta files download.**

Input can be entered manually or as CSV/TSV file:

* by choosing the manual input you’ll display a list of helpful TAGS to build your query

* otherwise if you choose to upload a CSV/TSV file,  you must specify which kind of list this one contains (Organism name list, taxonomy ID list or accession list); you may also narrowing down the search supplementing your query with genes or other details.

__Note__: list values must be placed within a single column (the first one).

When organism name list, Ids or accessions are less than 2500 the search key is composed by a single query, otherwise query will be splitted in groups of 2500 generating temporary files, which would be deleted at the end of the process.

Output file (standard format) will be saved in Download folder.

Available formats:
* multi-FASTA file format (NCBI standard format, with header followed by nucleotide sequence, accession, name, code, gene)
* TSV format


**Module 2: Statistical module; scatter plot and world map from taxonomy files or queries**

Module2’s required data, can be uploaded in several ways:
* manually, by entering a query;
* uploading module1’s taxonomy output file;
* uploading file from external sources containing taxonomy lists.

It's possible to do 3 types of interactive plots:
* scatterplot: uses taxonomy as input to produce a graph that indicates the quantity of each individual taxonomic unit and which taxons are present;
![alt text](https://github.com/albertobrusati/extaxsi_alberto/blob/master/aves%20scatterplot%20COX1.png)
* sunburst: uses taxonomy but creates an expansion pie that allows to explore taxonomy in depth with less weight on the quantity of each individual taxonomic unit;
![alt text](https://github.com/albertobrusati/extaxsi_alberto/blob/master/sunburst%20odonata.png)
* worldmap: uses the phylogeographic data to produce a map indicating the position of each individual species;
![alt text](https://github.com/albertobrusati/extaxsi_alberto/blob/master/worldmap.png)  

Output file format: html

Output file folder: download

**Module3:  Taxonomy ID converter.**

You can convert taxonomy ids from a file or by manual input into full 6 ranks taxonomy.

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

Contact information for the distributor or programmer
-----------------------------------------------------
* Name: __Anna Sandionigi__
* e-mail: __anna.sandionigi@unimib.it__
* Name: __Giulia Agostinetto__
* e-mail: __g.agostinetto@campus.unimib.it__

Credits and acknowledgments
---------------------------

Throughout the creation of this tool, we relied heavily on contributions from professors and college students.  Their input was invaluable, and we want to take a moment to thank them and recognize them for all of their hard work:
* __Dario Pescini__, Associate professor of University of Milan Bicocca.
* __Giulia Agostinetto__, PhD student of University of Milan Bicocca.
* __Adam Chahed__, student of University of Milan Bicocca.
* __Alberto Brusati__, student of University of Milan Bicocca.
* __Elena Parladori__, student of university of Milan Bicocca.
