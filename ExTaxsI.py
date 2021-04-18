#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ExTaxsI.py
#
#  Copyright 2019 Adam <a.chahed@campus.unimib.it>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

import time
import datetime
import os
import sys
import re
import http.client as http
import io

# To parse and use setting file
from configparser import ConfigParser

# For log files
import logging

# Needed to elaborate csv files
import csv
import glob
import pandas as pd

# This module implements container datatypes providing alternatives to
# Python’s general purpose built-in containers as counter
import collections
from collections import OrderedDict
from collections import Counter

# Biopython
from Bio import Entrez
from Bio.Entrez import Parser as ps
# from Bio import SeqIO

# import 1plotly.plotly as py
import plotly.graph_objects as go
from plotly.offline import plot

# Needed to search through xml files
import xml.etree.cElementTree as ET

#
from urllib.error import HTTPError

# To create a local taxonomy library to parse downloaded accessions
from ete3 import NCBITaxa
import sqlite3

# To check which OS we are using
import platform

# To give tab functionality while searching for files in bash (Not working in windows)
if platform.system() != "Windows":
    import readline

import numpy as np

import ETI_Lib as ETIL
import ETI_Lib_DB as ETIDB
import ETI_Lib_Interface as ETILI
import ETI_Lib_Statistics as ETILSTAT

ETIL.logInit()

'''
[ACCN]    Accession
[ALL]     All Fields
[AUTH]    Author
[GPRJ]    BioProject
[ECNO]    EC/RN Number
[FKEY]    Feature key
[FILT]    Filter
[GENE]    Gene Name
[JOUR]    Journal
[KYWD]    Keyword
[MLWT]    Molecular Weight
[ORGN]    Organism
[PACC]    Primary Accession
[PROP]    Properties
[PROT]    Protein Name
[SQID]    SeqID String
[SLEN]    Sequence Length
[SUBS]    Substance Name
[WORD]    Text Word
[TITL]    Title
[UID]     UID
'''



# # MAIN # #
# readline doesn't exist in windows so we don't use it in case we are in windows
if platform.system() != "Windows":
    readline.parse_and_bind('tab: complete')
    readline.set_completer_delims(' \t\n')

# Initializing
config = ConfigParser()

# read values from config file
config.read('settings.ini')

# Here goes your email so NCBI can contact you in case of necessity
logging.info(" Getting entrez_mail parameter from config file")
Entrez.email = config.get('parameters', 'entrez_email')

if Entrez.email == "example@example.com":
    print("Please enter your email in the settings.ini file")
    logging.info(" entrez_mail parameter is not defined (default value)")

else:
    logging.info(" Getting entrez_mail parameter completed")

# directory is the path where your downloads will go
logging.info(" Getting download_path parameter from config file")
directory = config.get('parameters', 'download_path')
logging.info(" Getting download_path parameter completed")

# Here goes your api key so
logging.info(" Getting api_key parameter from config file")
if config.get('parameters', 'api_key') == "none":
    logging.info(
        " No api_key parameter entered, for best results enter your ncbi api key in the settings.ini file"
    )

else:
    Entrez.api_key = config.get('parameters', 'api_key')
    logging.info(" Getting api_key parameter completed")

# Function to create the chosen download path
if not ETIL.create_folder(directory):
    directory = "./"

# Initialize local taxonomy database and, in case it isn't found, downloading and creating a local database
home = os.path.expanduser('~')

try:
    taxa_database_stat = os.stat(
        home +
        '/.etetoolkit/taxa.sqlite')  # Checking last update of the database
    taxa_database_time = taxa_database_stat.st_mtime
    print("\nTaxonomy database last update is {}".format(
        datetime.datetime.fromtimestamp(taxa_database_time)))
    taxa_database_update = str(
        input(
            "Do you want to update NCBI Taxonomy database? (Y > yes or N > no) "
        ))
    if taxa_database_update in ("Y", "y", "yes", "Yes", "YES"):
        print("Updating your local taxonomy database....")
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        print("\nUpdate completed!")
        time.sleep(2)

except FileNotFoundError:
    pass

ncbi = NCBITaxa()

while True:
    try:
        ETILI.clear()
        ETILI.main_menu()
    except KeyboardInterrupt:
        print(ETILI.color.RED + ETILI.color.BOLD +
              "\n\n<<KEYBOARD INTERRUPT HAS BEEN CAUGHT>>" + ETILI.color.END)
        print(
            "\nDo you want to go back to the main menu or you want to close the tool?"
        )
        if input(">> Enter 1 for main menu \n Enter anything else to Exit: "
                 ) == "1":
            continue
        else:
            print("\nClosing and saving all the files...")
            print("See you again!!")
            break
