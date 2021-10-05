#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ExTaxsI library.

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
import numpy as np
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

def load_configurations(entrez_email,api_key,taxa_database_update='no'):
    global directory
    global ncbi

    # To give tab functionality while searching for files in bash (Not working in windows)
    if platform.system() != "Windows":
        import readline

    # Initialize log_folder variable
    try:
        # In case the log folder in the tool folder doesn't exist
        if not os.path.exists("./Log/"):
            # It will try to create it
            os.makedirs("./Log/")
            # In case everything goes well, log file will be created in the log folder
            # The name of the file will be the date and the time at which the tool has started
        log_folder = './{0}/{1}.{2}'.format("Log", str(datetime.datetime.now().strftime("%d-%m-%y_at_%H-%M")), "log")

    except OSError:
        # In case of denied permissions or other cases where we can't create the log folder, the file will be in the same
        # folder of the tool
        # The name of the file will be the date and the time at which the tool has started
        log_folder = '{0}.{1}'.format(str(datetime.datetime.now().strftime("%d-%m-%y_at_%H-%M")), "log")

    # Creating the log file
    logging.basicConfig(filename=log_folder, format='%(asctime)s  %(levelname)s -> %(message)s', level=logging.DEBUG)

    # # MAIN # #
    # readline doesn't exist in windows so we don't use it in case we are in windows
    if platform.system() != "Windows":
        readline.parse_and_bind('tab: complete')
        readline.set_completer_delims(' \t\n')

    Entrez.email = entrez_email

    # directory is the path where your downloads will go
    directory = './download/'

    # Here goes your api key so

    Entrez.api_key = api_key

    # Function to create the chosen download path
    if not create_folder(directory):
        directory = "./"

    # Initialize local taxonomy database and, in case it isn't found, downloading and creating a local database
    home = os.path.expanduser('~')

    try:
        taxa_database_stat = os.stat(home + '/.etetoolkit/taxa.sqlite')  # Checking last update of the database
        taxa_database_time = taxa_database_stat.st_mtime
        print("\nTaxonomy database last update is {}".format(datetime.datetime.fromtimestamp(taxa_database_time)))
        if taxa_database_update in ("Y", "y", "yes", "Yes", "YES"):
            print("Updating your local taxonomy database....")
            ncbi = NCBITaxa()
            ncbi.update_taxonomy_database()
            print("\nUpdate completed!")
            time.sleep(2)

    except FileNotFoundError:
        pass

    ncbi = NCBITaxa()

class color:
    # Colors
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'

    # Text format
    BOLD = '\033[1m'
    ITALIC = '\33[3m'
    UNDERLINE = '\033[4m'
    BLINK = '\33[5m'

    # Background colors
    BLACKBG = '\33[40m'
    REDBG = '\33[41m'
    GREENBG = '\33[42m'
    YELLOWBG = '\33[43m'
    BLUEBG = '\33[44m'
    VIOLETBG = '\33[45m'
    BEIGEBG = '\33[46m'
    WHITEBG = '\33[47m'

    # End of format\colors
    END = '\033[0m'

# Function for bar progress animation
def update_progress(work_needed, work_done):
    # scritto da adam non funziona:
    #if not is_integer(work_done) or work_done <= 0:
    #    return

    bar_len = 40  # Modify this to change the length of the progress bar
    status = ""  # Initializing the string that will give us the status of the progression
    progress = work_needed / work_done

    if isinstance(progress, int):
        progress = float(progress)

    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"

    if progress < 0:  # progress can't be negative
        progress = 0
        status = "Halt...\r\n"

    if progress >= 1:
        progress = 1
        status = "Done...\r\n"

    block = int(round(bar_len * progress))  # % value of the progress
    text = '\r{6}Loading:{7} {5}[{0}]{7} {1}% || {2}/{3} {6}{4}{7}'.format('█' * block + '~' * (bar_len - block),
                                                                           round(progress * 100, 1), work_needed,
                                                                           work_done, status, color.GREEN, color.BOLD,
                                                                           color.END)
    sys.stdout.write(text)
    sys.stdout.flush()

def clear():
    if platform.system() == "Windows":
        os.system('cls')
    else:
        os.system('tput reset')

def merge_gene_top10(search_term,path,parent,manual_input=False,file_input=False):
    ###this part will create a unique dataframe from the multiple files retrieved by list-file search:
    if file_input:
        final_df = pd.DataFrame()
        os.chdir(path)
        for file in glob.glob(str('*_gene_list.tsv')):
            df = pd.read_csv(file, sep='\t', header=None)  # import as dataFrame
            df_t = df.T  # calculate transpose
            df_t.columns = df_t.iloc[0]
            df_t = df_t.drop([0])
            taxid = file.split('_')[0]  # extract taxid
            df_t.insert(0, 'Taxid', taxid)  # add taxid column at 0 position
            final_df = final_df.append(df_t, ignore_index = True, sort=False) # add row in final_df
            final_df = final_df.fillna(0)
        final_df.to_csv('final_gene_df.tsv', index=False, sep='\t')

        #final df to plot:
        final_df = final_df.iloc[:, 2:]
        df_t = final_df.T
        df_t['count'] = df_t.sum(axis=1)
        if len(df_t) >= 10:
            df_t = df_t.head(10).sort_values(by=['count'], ascending= True)
            df_t = df_t.rename_axis(['index']) #error file generation if name of axis is =0
            df_t.reset_index(level=0, inplace=True)
            df_top10 = df_t[['index','count']]
            graph = top10_graph(' list ',df_top10['count'],df_top10['index'])
            graph.write_image('top10_genes.png')
        else:
            df_t = df_t.sort_values(by=['count'], ascending= True)
            df_t = df_t.rename_axis(['index'])
            df_t.reset_index(level=0, inplace=True)
            df_top = df_t[['index','count']]
            graph = top10_graph(' list ',df_top['count'],df_top['index'])
            graph.write_image('top_used_genes.png')
        print('\nfinal_gene_df.tsv and top-plot created\n\n')
        os.chdir(parent)

    elif manual_input:
        os.chdir(path)
        dfplot = pd.read_csv(str(search_term.replace("(", "").replace(")", ""))+'_gene_list.tsv',  sep='\t', names=['gene','count'],skiprows=1)
        if len(dfplot) >= 10:
            df_top10 = dfplot.head(10).sort_values(by=['count'])
            graph = top10_graph(search_term,df_top10['count'],df_top10['gene'])
            graph.write_image('top10_genes_'+str(search_term.replace("(", "").replace(")", ""))+'.png')
        else:
            df_top10 = dfplot.sort_values(by=['count'])
            graph = top10_graph(search_term,df_top10['count'],df_top10['gene'])
            graph.write_image('top_used_genes_'+str(search_term.replace("(", "").replace(")", ""))+'.png')
        print('\n\nTOP-10 plotted!!!')
        os.chdir(parent)

def top10_graph(title, x_values, y_values):
    #y
    yval = y_values
    #x
    xval = x_values

    #plot
    fig=go.Figure([go.Bar(x=xval, y=yval,
                          text=xval,
                          textposition='auto',orientation='h')])

    #customize aspect
    fig.update_traces(overwrite=True,
                    marker={"color": xval,
                            "colorscale": 'blugrn'},
                                )

    fig.update_layout(title={'text':"Top10 "+title+": ",
                                    'y':0.95,
                                    'x':0.5,
                                    'xanchor': 'center',
                                    'yanchor': 'top'},
                                    autosize=False,
                                    width=900,
                                    height=950,
                                  #yaxis_title="",
                                  font=dict(family="Bodoni 72 Smallcaps",
                                            size=17,
                                            color='black'),
                                  paper_bgcolor='#ffffff',
                                  plot_bgcolor='#e6f0f5',
                                  bargap=0.4
                                )

    return fig

def download_gene_markers(counter, webenv, query_key, search_term, path, file):

    if counter == 0:
        return

    batch_size = 200
    genes_raw = []

    if counter > 400000:
        #print(webenv[0])
        indice = 0
        div_n = counter / 200000
        list_counter = []
        if div_n.is_integer():
            for num in range(0, int(div_n)):
                list_counter.append(200000)
        else:
            for num in range(0, int(div_n)):
                list_counter.append(200000)
            list_counter.append(counter-200000*int(div_n))

        #print(list_counter)

        for idx, mini_counter in enumerate(list_counter):

            for start in range(0, mini_counter, batch_size):
                logging.info(' Downloading data for gene markers: %i of %i' % (start + 1, mini_counter))
                data = efetch_call(start, "nuccore", mini_counter, webenv[indice], query_key, "gbc")

                if data is None:
                    print("\nMissing part!")
                    continue

                else:
                    try:
                        root = ET.fromstring(data)
                        INSDSeqs = root.findall('.//INSDSeq')

                        for INSDSeq in INSDSeqs:  # Iter every entry of our download
                            INSDQualifiers = INSDSeq.findall(
                                './/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier')

                            no_double_gene = []

                            for INSDQualifier in INSDQualifiers:
                                '''
                                if INSDQualifier[0].text == 'organism':
                                o anche accession per sapere da dove son stati presi i dati
                                    organism = str(INSDQualifier[1].text)
                                '''

                                if INSDQualifier[0].text == 'gene':

                                    if INSDQualifier[1].text in no_double_gene:
                                        continue
                                    else:
                                        no_double_gene.append(INSDQualifier[1].text)
                                        genes_raw.append(INSDQualifier[1].text)

                    except ET.ParseError as PE:
                        print("Error while parsing at %i \n" % start)
                        print("Error: %s" % PE)
                        root = None
                        INSDSeqs = None
                        INSDQualifiers = None
                        continue

                    except MemoryError as ME:
                        print("Error while parsing at %i \n" % start)
                        print("Error: %s" % ME)
                        root = None
                        INSDSeqs = None
                        INSDQualifiers = None
                        continue
            indice +=1

    else:
        for start in range(0, counter, batch_size):

            logging.info(' Downloading data for gene markers: %i of %i' % (start + 1, counter))
            data = efetch_call(start, "nuccore", counter, webenv, query_key, "gbc")

            if data is None:
                print("\nMissing part!")
                continue

            else:
                try:
                    root = ET.fromstring(data)
                    INSDSeqs = root.findall('.//INSDSeq')

                    for INSDSeq in INSDSeqs:  # Iter every entry of our download
                        INSDQualifiers = INSDSeq.findall(
                            './/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier')

                        no_double_gene = []

                        for INSDQualifier in INSDQualifiers:
                            '''
                            if INSDQualifier[0].text == 'organism':
                            o anche accession per sapere da dove son stati presi i dati
                                organism = str(INSDQualifier[1].text)
                            '''

                            if INSDQualifier[0].text == 'gene':

                                if INSDQualifier[1].text in no_double_gene:
                                    continue
                                else:
                                    no_double_gene.append(INSDQualifier[1].text)
                                    genes_raw.append(INSDQualifier[1].text)

                except ET.ParseError as PE:
                    print("Error while parsing at %i \n" % start)
                    print("Error: %s" % PE)
                    root = None
                    INSDSeqs = None
                    INSDQualifiers = None
                    continue

                except MemoryError as ME:
                    print("Error while parsing at %i \n" % start)
                    print("Error: %s" % ME)
                    root = None
                    INSDSeqs = None
                    INSDQualifiers = None
                    continue

    genes = []

    #  https://www.genecards.org/cgi-bin/carddisp.pl?gene=MT-CYB
    CYTB = ["Mitochondrially Encoded Cytochrome B", "Cytochrome B", "MTCYB", "CYTB",
            "Ubiquinol-Cytochrome-C Reductase Complex Cytochrome B", "Cytochrome B-C1", "COB"]
    COI = ["COXI", "Cytochrome C Oxidase I", "Cytochrome C Oxidase Subunit I", "Cytochrome C Oxidase Subunit 1", "COI",
            "COX1", "MTCO1", "Mitochondrially Encoded Cytochrome C Oxidase I", "Cytochrome C Oxidase Subunit I",
           "CO1"]
    COII = ["COXII", "Cytochrome C Oxidase II", "Cytochrome C Oxidase Subunit II", "Cytochrome C Oxidase Subunit 2",
            "COX2", "MTCO2", "Mitochondrially Encoded Cytochrome C Oxidase II","Cytochrome C Oxidase Subunit II",
            "COII", "CO2"]
    COIII = ["COXIII", "Cytochrome C Oxidase III", "Cytochrome C Oxidase Subunit III", "COX3", "MTCO3", "COIII",
             "Cytochrome C Oxidase Subunit 3", "Mitochondrially Encoded Cytochrome C Oxidase III", "CO3",
             "Cytochrome C Oxidase Subunit III"]
    ND1 = ["MTND1", "NADH1", "NADH Dehydrogenase 1", "Complex I ND1 Subunit", "NADH Dehydrogenase Subunit 1", "ND1",
           "NADH Ubiquinone Oxidoreductase Chain 1", "Mitochondrially Encoded NADH Dehydrogenase 1", "NAD1"]
    ND2 = ["MTND2", "NADH2", "NADH Dehydrogenase 2", "Complex I ND2 Subunit", "NADH Dehydrogenase Subunit 2", "ND2",
           "NADH Ubiquinone Oxidoreductase Chain 2", "Mitochondrially Encoded NADH Dehydrogenase 2", "NAD2"]
    ND3 = ["MTND3", "NADH3", "NADH Dehydrogenase 3", "Complex I ND3 Subunit", "NADH Dehydrogenase Subunit 3", "ND3",
           "NADH Ubiquinone Oxidoreductase Chain 3", "Mitochondrially Encoded NADH Dehydrogenase 3", "NAD3"]
    ND4 = ["MTND4", "NADH4", "NADH Dehydrogenase 4", "Complex I ND4 Subunit", "NADH Dehydrogenase Subunit 4", "ND4",
           "NADH Ubiquinone Oxidoreductase Chain 4", "Mitochondrially Encoded NADH Dehydrogenase 4", "NAD4"]
    ND5 = ["MTND5", "NADH5", "NADH Dehydrogenase 5", "Complex I ND5 Subunit", "NADH Dehydrogenase Subunit 5", "ND5",
           "NADH Ubiquinone Oxidoreductase Chain 5", "Mitochondrially Encoded NADH Dehydrogenase 5", "NAD5"]
    ND6 = ["MTND6", "NADH6", "NADH Dehydrogenase 6", "Complex I ND6 Subunit", "NADH Dehydrogenase Subunit 6", "ND6",
           "NADH Ubiquinone Oxidoreductase Chain 6", "Mitochondrially Encoded NADH Dehydrogenase 6", "NAD6"]
    ATP6 = ["ATP6", "MTATP6", "Mitochondrially Encoded ATP Synthase Membrane Subunit 6", "ATP Synthase 6", "ATPase6",
            "ATP Synthase Subunit A", "Mitochondrially Encoded ATP Synthase 6", "ATP Synthase F0 Subunit 6"]
    ATP8 = ["ATP8", "MTATP8", "Mitochondrially Encoded ATP Synthase Membrane Subunit A6L", "ATP Synthase 8", "ATPase8",
            "ATP Synthase Subunit A", "Mitochondrially Encoded ATP Synthase 8", "ATP Synthase F0 Subunit 8", "A6L"]

    for gene_raw in genes_raw:
        gene = re.sub(r'[^a-zA-Z0-9]', "", gene_raw.upper())

        if gene in CYTB:
            genes.append("CYTB")
            continue

        elif gene in ATP6:
            genes.append("ATP6")
            continue

        elif gene in ATP8:
            genes.append("ATP8")
            continue

        elif gene in COIII:
            genes.append("COIII")
            continue

        elif gene in COII:
            genes.append("COII")
            continue

        elif gene in COI:
            genes.append("COI")
            continue

        elif gene in ND1:
            genes.append("ND1")
            continue

        elif gene in ND2:
            genes.append("ND2")
            continue

        elif gene in ND3:
            genes.append("ND3")
            continue

        elif gene in ND4:
            genes.append("ND4")
            continue

        elif gene in ND5:
            genes.append("ND5")
            continue

        elif gene in ND6:
            genes.append("ND6")
            continue

        elif "16S" in gene:
            genes.append("16S")
            continue

        elif "18S" in gene:
            genes.append("18S")
            continue

        else:
            genes.append(gene_raw.upper())

    gene_array = np.array(genes)
    gene_counter_list = collections.Counter(gene_array)
    update_progress(counter, counter)

    if file is None:
        file = "{0}{1}_gene_list.tsv".format(path, search_term.replace("(", "").replace(")", ""))

    out_handle = open(file, "w")
    out_handle.write("n_records_found\t{0}\n".format(counter))

    for gene in OrderedDict(gene_counter_list.most_common()):
        out_handle.write("{0}\t{1}\n".format(gene, gene_counter_list[gene]))

    out_handle.close()
    print("File created at %s\n" % file)

# Creates the path where the files will be downloaded
def create_folder(folder_path):
    logging.info(' Searching or creating the download path')

    try:
        if not os.path.exists(folder_path):  # If the folder doesn't exist it will make it
            os.makedirs(folder_path)
            logging.info(' Path created!')
            return True

        else:
            logging.info(' Path already exist')
            return True
    except OSError:
        print('Error: Creating directory ' + folder_path, '\nFiles will be saved in the tool directory')
        logging.warning(
            ' Error at creating directory << %s >> try to change your path in the settings file' % folder_path)
        return False

# Rename files name
def rename_file(folder_path, rename_me, extension):

    while "/" in rename_me:

        if "/" in rename_me:
            rename_me = rename_me[rename_me.find("/") + 1:]

    if not os.path.exists(folder_path):
        rename_me = '{0}{1}'.format(rename_me.replace(" ", "_").replace("(", "").replace(")", ""), extension)

    else:
        rename_me = '{0}{1}{2}'.format(folder_path, rename_me.replace(" ", "_").replace("(", "").replace(")", ""),
                                       extension)

    return rename_me

def efetch_call(start, db, counter_id, webenv, query_key, rettype):

    batch_size = 200
    update_progress(start, counter_id)  # Download animation bar function

    server_errors = 0
    while server_errors < 11:  # More than one attempts every time there is a network

        try:

            fetch_handle = Entrez.efetch(db=db,
                                         retstart=start,
                                         retmax=batch_size,
                                         webenv=webenv,
                                         query_key=query_key,
                                         rettype=rettype)

            if rettype == "fasta" or rettype == "gbc":
                data = fetch_handle.read()

            else:
                data = Entrez.read(fetch_handle, validate=False)

            fetch_handle.close()
            break

        except HTTPError as err:
            server_errors += 1

            if 500 <= err.code <= 599:
                logging.warning(' Connection error to the server %s' % err)
                print(color.RED + '\nConnection error: %s' % err + color.END)
                logging.warning(" Attempt %i of 10" % server_errors)
                print("Attempt %i of 10" % server_errors)
                time.sleep(10)
                if server_errors > 6:
                    print("Maybe NCBI's servers are too busy, if it still doesn't work, try later")

            else:
                print("\nError in the connection: %s" % err)
                logging.error("Exception occurred", exc_info=True)
                print("Attempt %i of 10" % server_errors)
                time.sleep(20)

        except OSError as err:
            server_errors += 1
            print("\nError in the connection: %s" % err)
            logging.error("Exception occurred", exc_info=True)
            print("Attempt %i of 10" % server_errors)
            time.sleep(20)

        except http.IncompleteRead as err:
            print("Error: %s" % err)
            time.sleep(5)

        except RuntimeError as err:
            print("Error: %s" % err)
            print("Retrying in a few seconds")
            time.sleep(10)

    if server_errors >= 11:
        logging.error(" Couldn't fetch part of the data list at: %i" % start)
        print("\nMissing part")
        logging.error("Try again starting from %i" % start)
        fetch_handle.close()
        return None

    else:
        return data

# Download dei FASTA
def download_fasta(counter_id, webenv, query_key, query, folder_path, file_input):
    batch_size = 200  # Batch_size value limit our queries to NCBI so we don't get blacklisted
    if file_input is None:
        file_name = rename_file(folder_path, query, ".fasta")  # Formatting the file name
        out_handle = open(file_name, "w")  # Creating the fasta file, where we will save our downloads

    logging.info('Downloading Fasta files')
    missing_part = 0

    print("Downloading fasta sequence")
    if counter_id > 400000:
        #print(webenv[0])
        indice = 0
        div_n = counter_id / 200000
        list_counter = []
        if div_n.is_integer():
            for num in range(0, int(div_n)):
                list_counter.append(200000)
        else:
            for num in range(0, int(div_n)):
                list_counter.append(200000)
            list_counter.append(counter_id-200000*int(div_n))

        #print(list_counter)

        for idx, mini_counter in enumerate(list_counter):

            for start in range(0, mini_counter, batch_size):
                logging.info('Downloading fasta sequence: %i of %i' % (start + 1, mini_counter))
                data = efetch_call(start, "nucleotide", mini_counter, webenv[indice], query_key, "fasta")

                if data is None:
                    missing_part = 1

                elif file_input is None:
                    out_handle.write(data)

                else:
                    with open(file_input, 'a+') as a_writer:
                        a_writer.write("\n{0}".format(data))
            indice += 1
    else:
        for start in range(0, counter_id, batch_size):
            logging.info('Downloading fasta sequence: %i of %i' % (start + 1, counter_id))
            data = efetch_call(start, "nucleotide", counter_id, webenv, query_key, "fasta")

            if data is None:
                missing_part = 1

            elif file_input is None:
                out_handle.write(data)

            else:
                with open(file_input, 'a+') as a_writer:
                    a_writer.write("\n{0}".format(data))

    if file_input is None:
        out_handle.close()

    else:
        a_writer.close()

    update_progress(counter_id, counter_id)

    if missing_part == 0:
        logging.info(' Fasta file has been created')
        print("\n ---- Fasta file has been created. ----\n\n")

    else:
        print("\n ---- Fasta file has been created but some errors occurred. ----")
        print(" ---- Check log file for more details. ----\n\n")

# Create a file with ACCESSION iD and six level of Taxonomy (Phylum,Class,Order,Family,Genus,Species)
def download_accession_taxonomy(counter_id, webenv, query_key, query, folder_path, file_input):

    logging.info(' Creating accession file')
    if file_input is None:
        file_name = query
        file_name = rename_file(folder_path, file_name, "_taxonomy.tsv")
        out_handle = open(file_name, "w")

    batch_size = 200
    missing_part = 0

    print("Downloading accession with taxonomy file")
    if counter_id > 400000:
        #print(webenv[0])
        indice = 0
        div_n = counter_id / 200000
        list_counter = []
        if div_n.is_integer():
            for num in range(0, int(div_n)):
                list_counter.append(200000)
        else:
            for num in range(0, int(div_n)):
                list_counter.append(200000)
            list_counter.append(counter_id-200000*int(div_n))

        #print(list_counter)

        for idx, mini_counter in enumerate(list_counter):

            for start in range(0, mini_counter, batch_size):
                logging.info('\nDownloading accession with taxonomy file: %i of %i' % (start + 1, mini_counter))

                data = efetch_call(start, "nucleotide", mini_counter, webenv[indice], query_key, "gpc")


                if data is None:
                    missing_part = 1
                    continue

                else:

                    for n in range(len(data)):  # For every downloaded entry
                        tax_id = None

                        # if taxonomy id is not found, it will give NA for every taxonomy level
                        try:
                            tax_id = ncbi.get_name_translator([data[n]["INSDSeq_organism"]])

                            if tax_id:
                                tax_id = tax_id[data[n]["INSDSeq_organism"]][0]

                            else:
                                tax_id = None

                        except ValueError:
                            logging.warning(" There is no taxonomy for %s \n" +
                                            " searching through taxID" % data[n]["INSDSeq_organism"])

                        except sqlite3.OperationalError:
                            logging.warning(" There is no taxonomy for %s \n " +
                                            "searching through taxID" % data[n]["INSDSeq_organism"])

                        if tax_id is None:
                            for z in range(len(data[n]["INSDSeq_feature-table"][0]["INSDFeature_quals"])):
                                # Searching for taxonomy id in the xml type data
                                try:
                                    if "taxon" in data[n]["INSDSeq_feature-table"][0]["INSDFeature_quals"][z]["INSDQualifier_value"]:
                                        tax_id = data[n]["INSDSeq_feature-table"][0]["INSDFeature_quals"][z][
                                            "INSDQualifier_value"].replace("taxon:", "")
                                        break

                                except KeyError:
                                    continue

                        if tax_id is None:
                            logging.error(
                                " There is no taxonomy for %s, ncbi.get_lineage exception" % data[n]["INSDSeq_organism"],
                                exc_info=True)
                            taxa_list = "NA;NA;NA;NA;NA"
                            accession_list = "{0}\t{1};{2}{3}".format(data[n]["INSDSeq_accession-version"], taxa_list,
                                                                       data[n]["INSDSeq_organism"], "\n")
                            if file_input is None:
                                out_handle.write(accession_list)

                            else:
                                with open(file_input, 'a+') as a_writer:
                                    a_writer.write(accession_list)
                            continue

                        try:
                            lineage = ncbi.get_lineage(tax_id)  # Getting entire taxonomy from its taxonomy id

                        except ValueError:  # In case its id is not found
                            logging.error(
                                " There is no taxonomy for %s, ncbi.get_lineage exception" % data[n]["INSDSeq_organism"],
                                exc_info=True)
                            taxa_list = "NA;NA;NA;NA;NA"
                            accession_list = "{0}\t{1};{2}{3}".format(data[n]["INSDSeq_accession-version"], taxa_list,
                                                                       data[n]["INSDSeq_organism"], "\n")
                            if file_input is None:
                                out_handle.write(accession_list)

                            else:
                                with open(file_input, 'a+') as a_writer:
                                    a_writer.write(accession_list)
                            continue

                        except sqlite3.OperationalError as SQ:
                            print("sqlite3 error")
                            print(SQ)
                            print(str(tax_id) + "\n")
                            time.sleep(5)
                            continue

                        phylum = clas = order = family = genus = "NA"  # Initializing

                        if lineage is not None:
                            for z in range(len(lineage)):
                                lineage_rank = ncbi.get_rank([lineage[z]])

                                # Checking the rank and getting their name
                                if "phylum" == lineage_rank[lineage[z]]:
                                    rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                    phylum = rank_tmp[lineage[z]]

                                if "class" == lineage_rank[lineage[z]]:
                                    rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                    clas = rank_tmp[lineage[z]]

                                if "order" == lineage_rank[lineage[z]]:
                                    rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                    order = rank_tmp[lineage[z]]

                                if "family" == lineage_rank[lineage[z]]:
                                    rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                    family = rank_tmp[lineage[z]]

                                if "genus" == lineage_rank[lineage[z]]:
                                    rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                    genus = rank_tmp[lineage[z]]

                        taxa_list = ";".join([phylum, clas, order, family, genus])
                        accession_list = "{0}\t{1};{2}{3}".format(data[n]["INSDSeq_accession-version"], taxa_list,
                                                                  data[n]["INSDSeq_organism"], "\n")
                        if file_input is None:
                            out_handle.write(accession_list)

                        else:
                            with open(file_input, 'a+') as a_writer:
                                a_writer.write(accession_list)

            indice += 1

        if file_input is None:
            out_handle.close()

        else:
            a_writer.close()

        update_progress(counter_id, counter_id)
        logging.info(' Accession with taxonomy file has been created')

    else:
        for start in range(0, counter_id, batch_size):
            logging.info(' Downloading accession with taxonomy file: %i of %i' % (start + 1, counter_id))

            data = efetch_call(start, "nucleotide", counter_id, webenv, query_key, "gpc")

            if data is None:
                missing_part = 1
                continue

            else:

                for n in range(len(data)):  # For every downloaded entry
                    tax_id = None

                    # if taxonomy id is not found, it will give NA for every taxonomy level
                    try:
                        tax_id = ncbi.get_name_translator([data[n]["INSDSeq_organism"]])

                        if tax_id:
                            tax_id = tax_id[data[n]["INSDSeq_organism"]][0]

                        else:
                            tax_id = None

                    except ValueError:
                        logging.warning(" There is no taxonomy for %s \n" +
                                        " searching through taxID" % data[n]["INSDSeq_organism"])

                    except sqlite3.OperationalError:
                        logging.warning(" There is no taxonomy for %s \n " +
                                        "searching through taxID" % data[n]["INSDSeq_organism"])

                    if tax_id is None:
                        for z in range(len(data[n]["INSDSeq_feature-table"][0]["INSDFeature_quals"])):
                            # Searching for taxonomy id in the xml type data
                            try:
                                if "taxon" in data[n]["INSDSeq_feature-table"][0]["INSDFeature_quals"][z]["INSDQualifier_value"]:
                                    tax_id = data[n]["INSDSeq_feature-table"][0]["INSDFeature_quals"][z][
                                        "INSDQualifier_value"].replace("taxon:", "")
                                    break

                            except KeyError:
                                continue

                    if tax_id is None:
                        logging.error(
                            " There is no taxonomy for %s, ncbi.get_lineage exception" % data[n]["INSDSeq_organism"],
                            exc_info=True)
                        taxa_list = "NA;NA;NA;NA;NA"
                        accession_list = "{0}\t{1};{2}{3}".format(data[n]["INSDSeq_accession-version"], taxa_list,
                                                                   data[n]["INSDSeq_organism"], "\n")
                        if file_input is None:
                            out_handle.write(accession_list)

                        else:
                            with open(file_input, 'a+') as a_writer:
                                a_writer.write(accession_list)
                        continue

                    try:
                        lineage = ncbi.get_lineage(tax_id)  # Getting entire taxonomy from its taxonomy id

                    except ValueError:  # In case its id is not found
                        logging.error(
                            " There is no taxonomy for %s, ncbi.get_lineage exception" % data[n]["INSDSeq_organism"],
                            exc_info=True)
                        taxa_list = "NA;NA;NA;NA;NA"
                        accession_list = "{0}\t{1};{2}{3}".format(data[n]["INSDSeq_accession-version"], taxa_list,
                                                                   data[n]["INSDSeq_organism"], "\n")
                        if file_input is None:
                            out_handle.write(accession_list)

                        else:
                            with open(file_input, 'a+') as a_writer:
                                a_writer.write(accession_list)
                        continue

                    except sqlite3.OperationalError as SQ:
                        print("sqlite3 error")
                        print(SQ)
                        print(str(tax_id) + "\n")
                        time.sleep(5)
                        continue

                    phylum = clas = order = family = genus = "NA"  # Initializing

                    if lineage is not None:
                        for z in range(len(lineage)):
                            lineage_rank = ncbi.get_rank([lineage[z]])

                            # Checking the rank and getting their name
                            if "phylum" == lineage_rank[lineage[z]]:
                                rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                phylum = rank_tmp[lineage[z]]

                            if "class" == lineage_rank[lineage[z]]:
                                rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                clas = rank_tmp[lineage[z]]

                            if "order" == lineage_rank[lineage[z]]:
                                rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                order = rank_tmp[lineage[z]]

                            if "family" == lineage_rank[lineage[z]]:
                                rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                family = rank_tmp[lineage[z]]

                            if "genus" == lineage_rank[lineage[z]]:
                                rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                                genus = rank_tmp[lineage[z]]

                    taxa_list = ";".join([phylum, clas, order, family, genus])
                    accession_list = "{0}\t{1};{2}{3}".format(data[n]["INSDSeq_accession-version"], taxa_list,
                                                              data[n]["INSDSeq_organism"], "\n")
                    if file_input is None:
                        out_handle.write(accession_list)

                    else:
                        with open(file_input, 'a+') as a_writer:
                            a_writer.write(accession_list)

        if file_input is None:
            out_handle.close()

        else:
            a_writer.close()

        update_progress(counter_id, counter_id)
        logging.info(' Accession with taxonomy file has been created')

    if missing_part == 0:
        print("\n ---- Accession with taxonomy file has been created. ----\n")

    else:
        print("\n ---- Accession with taxonomy file has been created but some errors occurred. ----")
        print(" ---- Check log file for more details. ---- \n")
    if file_input is None:
        return file_name
    else:
        return

def db_creation(text_search = None,
                file_search = None,
                input_file_type = None,
                additional_query = 0,
                fasta_output = False,
                accession_taxonomy_output = False,
                marker_output = False,
                top10_plot = False,
                enrich_output = False):

    if text_search is not None:
        counter_query = 0
        while counter_query == 0:
            server_errors = 0

            while server_errors < 11:
                server_errors += 1

                try:
                    handle = Entrez.esearch(db="nucleotide", term=text_search, usehistory="y", idtype="acc")
                    break

                except HTTPError as err:
                    server_errors += 1

                    if 500 <= err.code <= 599:
                        logging.warning(' Connection error to the server %s' % err)
                        print('Connection error to the server %s' % err)
                        logging.warning(" Attempt %i of 10" % server_errors)
                        print("Attempt %i of 10" % server_errors)
                        time.sleep(10)
                        if server_errors > 6:
                            print("Maybe NCBI's servers are too busy, if it still doesn't work, try later")

                    else:
                        print("\nError in the connection: %s" % err)
                        logging.error("Exception occurred", exc_info=True)
                        print("Attempt %i of 10" % server_errors)
                        time.sleep(20)

                except OSError as err:
                    print("\nError in the connection: %s" % err)
                    logging.error("Exception occurred", exc_info=True)
                    print("Attempt %i of 10" % server_errors)
                    time.sleep(20)

                except http.IncompleteRead as err:
                    print("Error: %s" % err)
                    time.sleep(5)

                except RuntimeError as err:
                    print("Error: %s" % err)
                    print("Retrying in a few seconds")
                    time.sleep(10)

            if server_errors >= 11:
                print("Possible maintenance at NCBI servers if the problem is not in your connection, " +
                      "try again some minutes later")
                logging.warning("Connection impossible with NCBI, can't request %s" % text_search)
                continue

            results = Entrez.read(handle, validate=False)

            # results["count"] contains the results counts found with the query
            # This number is needed to define how many ID's we will search and download
            # Default value is 20 means NCBI will let us download only a bunch of data instead of the full list
            counter_query = int(results["Count"])

            if counter_query == 0:
                print("\nNo results found with this query: %s\nMaybe a typo?\n" % text_search)
                continue
            else:
                print("\nThe number of sequence found is: {0}{2}{1}".format(color.GREEN,
                                                                            color.END,
                                                                            counter_query))
                if counter_query > 400000:
                    print(color.RED + "High number of records!" + color.END +
                          " It can take a while to load..")

            server_errors = 0

            while server_errors < 8:
                server_errors += 1

                try:
                    handle = Entrez.esearch(db="nucleotide",
                                            term=text_search,
                                            retmax=counter_query)
                    results = Entrez.read(handle, validate=False)

                    id_list = results["IdList"]

                    if len(id_list) > 400000:
                        webenv = []
                        query_key = []

                        for index in range(0, len(id_list), 200000):
                            ###
                            print('Splitting query: ')
                            print(index)
                            ###
                            mini_id_list = []
                            end = min(index + 200000, len(id_list))

                            for id in range(index, end):
                                mini_id_list.append(id_list[id])
                            #print(len(mini_id_list))
                            # Using epost we will lighten the burden on Entrez
                            post_xml = Entrez.epost("nucleotide",
                                                    id=",".join(mini_id_list))
                            results = Entrez.read(post_xml, validate=False)

                            # Variable needed to use ncbi history feature
                            webenv.append(results["WebEnv"])
                            #print(webenv)
                            query_key.append(results["QueryKey"])
                        #print(webenv)
                        #print(query_key)
                    else:
                        post_xml = Entrez.epost("nucleotide",
                                                id=",".join(id_list))
                        results = Entrez.read(post_xml, validate=False)
                        # Variable needed to use ncbi history feature
                        webenv = results["WebEnv"]
                        query_key = results["QueryKey"]

                    break

                except HTTPError as err:

                    if 500 <= err.code <= 599:
                        logging.warning(' Connection error to the server %s' % err)
                        print('\nConnection error to the server %s' % err)
                        logging.warning(" Attempt %i of 7" % server_errors)
                        print("Attempt %i of 7" % server_errors)
                        time.sleep(10)
                        if server_errors > 6:
                            print("Maybe NCBI's servers are too busy, if it still doesn't work, try later")

                    else:
                        print("\nError in the connection: %s" % err)
                        logging.error("Exception occurred", exc_info=True)
                        print("Attempt %i of 7" % server_errors)
                        time.sleep(20)

                except OSError as err:
                    print("\nError in the connection: %s" % err)
                    logging.error("Exception occurred", exc_info=True)
                    print("Attempt %i of 7" % server_errors)
                    time.sleep(20)

                except http.IncompleteRead as err:
                    print("Error: %s" % err)
                    time.sleep(5)

                except RuntimeError as err:
                    print("Error: %s" % err)
                    print("Retrying in a few seconds")
                    time.sleep(10)

                except ps.CorruptedXMLError as err:
                    print("XML error: %s" % err)

            if server_errors >= 8:
                print("Connection impossible with NCBI, can't request %s" % text_search)
                logging.warning("Connection impossible with NCBI, can't request %s" % text_search)
                continue
        if fasta_output:
            download_fasta(counter_query,
                           webenv,
                           query_key,
                           text_search,
                           directory,
                           None)
        if accession_taxonomy_output:
            download_accession_taxonomy(counter_query,
                                        webenv,
                                        query_key,
                                        text_search,
                                        directory,
                                        None)
        if marker_output:
            download_gene_markers(counter_query,
                                  webenv,
                                  query_key,
                                  text_search,
                                  directory,
                                  None,)
            parentdir = os.getcwd()
            merge_gene_top10(text_search,directory,parentdir, manual_input = top10_plot)

        if enrich_output:
            download_enrich_file(counter_query,
                                 webenv,
                                 query_key,
                                 text_search,
                                 directory)

    if file_search is not None:

        if input_file_type in ("T", "t", "O", "o"):
            # POSSIBLE INPUT -> GENE IN THE COLUMN 0
            file_list = additional_query

        elif input_file_type in ("A", "a"):
            file_list = ["accession"]

        else:
            print("\n Parameter <input_file_type> must be defined for a file search!\nPlease choose one of the following: 'T', 'O', 'A' \n")
            quit()

        while True:
            csv_file_module_one = file_search
            print("\nLoading files...\n")

            try:
                logging.info("Loading csv file...")
                with open(csv_file_module_one, 'r') as r:

                    while csv_file_module_one[-3:] not in ('tsv', 'csv'):
                        print("File extension not CSV or TSV, please use the right format.")
                        csv_file_module_one = input("Enter your csv or tsv file path: ")

                    if csv_file_module_one[-3:] == 'csv':
                        file_csv_query = csv.reader(r, delimiter=',')

                    if csv_file_module_one[-3:] == 'tsv':
                        file_csv_query = csv.reader(r, delimiter='\t')

                    for line in file_csv_query:
                        file_list.append(line[0])

                    logging.info("Loading completed!")
                    break

            except IOError as err:
                logging.warning('No file found or permission denied, please check your file or location: %s' % err)
                print("No file found or permission denied, please check your file or location")
                print("File location: ", csv_file_module_one, "\n")

        while "/" in csv_file_module_one:
            csv_file_module_one = csv_file_module_one[csv_file_module_one.index("/") + 1 : ]


        if marker_output:
            first_element = True

            if create_folder("./download/{0}/".format(csv_file_module_one)):
                marker_folder = "./download/{0}/".format(csv_file_module_one)

            else:
                marker_folder = directory

            for element in file_list:

                if first_element:
                    first_element = False
                    continue

                if input_file_type in ("o", "O"):  # File list of Organism

                    if file_list[0] == "0":
                        search_index = element

                    else:
                        search_index = "({0} AND {1})".format(element, file_list[0])

                elif input_file_type in ("t", "T"):  # File list of TaxIDs

                    if file_list[0] == "0":
                        search_index = "{1}{0}{2}".format(element, "txid", "[Organism:noexp]")

                    else:
                        search_index = "({2} AND {1}{0}{3})".format(element,
                                                                    "txid",
                                                                    file_list[0],
                                                                    "[Organism:noexp]")

                else:  # File list of ACCESSIONS
                    search_index = [str(element)]

                # variable to prevent the stop of esearch if you get a momentary internet problem
                server_error = 0

                while server_error < 11:
                    server_error += 1

                    try:
                        handles = Entrez.esearch(db="nucleotide",
                                                 term=search_index,
                                                 usehistory="y",
                                                 idtype="acc")

                        result = Entrez.read(handles, validate=False)

                        # With counter_query we know how many IDs we found
                        counter_queries = int(result["Count"])

                        if counter_queries == 0:
                            print("\nNo results found with: \n",
                                  search_index, "\n")
                            logging.warning("No results found with: %s" % search_index)
                            break

                        # results["count"] contains the results counts found with the query
                        # This number is needed to define how many ID's we will search and download
                        # Default value is 20 and NCBI will let us download only 20 ids instead of the full list
                        handles = Entrez.esearch(db="nucleotide",
                                                 term=search_index,
                                                 retmax=counter_queries)

                        result = Entrez.read(handles, validate=False)
                        id_list = result["IdList"]
                        epost_xml = Entrez.epost("nucleotide", id=",".join(id_list))
                        result = Entrez.read(epost_xml, validate=False)
                        web_env = result["WebEnv"]
                        key = result["QueryKey"]
                        break

                    except HTTPError as err:
                        server_error += 1

                        if 500 <= err.code <= 599:
                            logging.warning(' Connection error to the server %s' % err)
                            print('\nConnection error to the server %s' % err)
                            print(search_index)
                            time.sleep(60)
                            logging.warning(" Attempt %i of 10" % server_error)
                            print("Attempt %i of 10" % server_error)
                            time.sleep(10)
                            if server_error > 6:
                                print("Maybe NCBI's servers are too busy, if it still doesn't work, try later")

                        else:
                            print("\nError in the connection: %s" % err)
                            logging.error("Exception occurred", exc_info=True)
                            print("Attempt %i of 10" % server_error)
                            time.sleep(20)

                    except OSError as err:
                        print("\nError in the connection: %s" % err)
                        logging.error("Exception occurred", exc_info=True)
                        print("Attempt %i of 10" % server_error)
                        time.sleep(20)

                    except http.IncompleteRead as err:
                        print("Error: %s" % err)
                        time.sleep(5)

                    except RuntimeError as err:
                        print("Error: %s" % err)
                        print("Retrying in a few seconds")
                        time.sleep(10)

                if server_error >= 11:
                    print("Connection impossible with NCBI, can't request %s" % search_index)
                    logging.warning("Connection impossible with NCBI, can't request %s" % search_index)
                    #Bloccare?

                file_name = rename_file(marker_folder,
                                        element,
                                        "_gene_list.tsv")
                print(file_name)

                if os.path.exists(file_name):
                    x = 0
                    overwrite = str(input(color.YELLOW +
                                          "Exist already a file with the same name, do you want to overwrite it? " +
                                          "(0 = yes, 1 = no) " +
                                          color.END))

                else:
                    overwrite = 1

                while os.path.exists(file_name) and overwrite != "0":
                    x += 1
                    file_name = rename_file(directory, element, "_gene_list_%i.tsv" % x)

                if counter_queries != 0:
                    download_gene_markers(counter_queries, web_env, key, search_index, directory, file_name)
                    time.sleep(2)

            parentdir = os.getcwd()
            merge_gene_top10(None, marker_folder,parentdir,file_input=top10_plot)

        if len(file_list) > 2500:  # NCBI rejects more than 2500 - 2700 terms in a unique call so we will split them
            search_list = [file_list[0]]
            for pos in range(1, 2500):
                search_list.append(file_list[pos])
        else:
            search_list = file_list

        if input_file_type in ("o", "O"):  # File list of Organism
            if search_list[0] == "0":
                search_index = "{0} OR {1}".format(search_list[1], search_list[2])
                for i in range(3, len(search_list)):
                    search_index = "{0} OR {1}".format(search_index, search_list[i])

            else:
                search_index = "({0} AND {1})".format(search_list[1], search_list[0])

                for i in range(2, len(search_list)):
                    search_index = "{2} OR ({0} AND {1})".format(search_list[i], search_list[0], search_index)

        elif input_file_type in ("t", "T"):  # File list of TaxIDs

            if search_list[0] == "0":
                search_index = "{2}{0}{3} OR {2}{1}{3}".format(search_list[1], search_list[2], "txid","[Organism:noexp]")
                for i in range(3, len(search_list)):
                    search_index = "{0} OR {2}{1}{3}".format(search_index, search_list[i], "txid", "[Organism:noexp]")

            else:
                search_index = "({2} AND {1}{0}{3})".format(search_list[1], "txid", search_list[0], "[Organism:noexp]")
                for i in range(2, len(search_list)):
                    search_index = "{0} OR ({3} AND {2}{1}{4})".format(search_index, search_list[i], "txid",
                                                                    search_list[0], "[Organism:noexp]")

        else:  # File list of ACCESSIONS
            search_index = [str(search_list[1])]

            for i in range(2, len(search_list)):
                search_index = "{0} or {1}".format(search_index, search_list[i])

        # variable to prevent the stop of esearch if you get a momentary internet problem
        server_error = 0

        while server_error < 11:

            try:
                handles = Entrez.esearch(db="nucleotide", term=search_index, usehistory="y", idtype="acc")
                break

            except HTTPError as err:
                server_error += 1

                if 500 <= err.code <= 599:
                    logging.warning(' Connection error to the server %s' % err)
                    print('\nConnection error to the server %s' % err)
                    print(search_index)
                    time.sleep(60)
                    logging.warning(" Attempt %i of 10" % server_error)
                    print("Attempt %i of 10" % server_error)
                    time.sleep(10)
                    if server_error > 6:
                        print("Maybe NCBI's servers are too busy, if it still doesn't work, try later")

                else:
                    print("\nError in the connection: %s" % err)
                    logging.error("Exception occurred", exc_info=True)
                    print("Attempt %i of 10" % server_error)
                    time.sleep(20)

            except OSError as err:
                print("\nError in the connection: %s" % err)
                logging.error("Exception occurred", exc_info=True)
                print("Attempt %i of 10" % server_error)
                time.sleep(20)

            except http.IncompleteRead as err:
                print(err)
                print("\n Retrying in a few seconds")
                time.sleep(10)

            except RuntimeError as err:
                print(err)
                print("\n Retrying in a few seconds")
                time.sleep(10)

        if server_error >= 11:
            print("Connection impossible with NCBI, can't request %s" % search_index)
            logging.warning("Connection impossible with NCBI, can't request %s" % search_index)
            quit()

        result = Entrez.read(handles, validate=False)
        counter_queries = int(result["Count"])  # With counter_query we know how many IDs we found

        if counter_queries == 0:
            print("\nNo results found with: \n", search_index, "\n")
            logging.warning("No results found with: %s" % search_index)

        # results["count"] contains the results counts found with the query
        # This number is needed to define how many ID's we will search and download
        # Default value is 20 so without it NCBI will let us download only a bunch of data instead of the full list
        handles = Entrez.esearch(db="nucleotide", term=search_index, retmax=counter_queries)

        result = Entrez.read(handles, validate=False)
        id_list = result["IdList"]
        epost_xml = Entrez.epost("nucleotide", id=",".join(id_list))
        result = Entrez.read(epost_xml, validate=False)
        web_env = result["WebEnv"]
        key = result["QueryKey"]

        if fasta_output:
            file_name = rename_file(directory, csv_file_module_one, ".fasta")

            if os.path.exists(file_name):
                x = 0
                overwrite = str(input(color.YELLOW +
                                  "0 => yes, 1 => no" +
                                  "Exist already a file with the same name, do you want to overwrite it? "
                                  + color.END))
            else:
                overwrite = 1

            while os.path.exists(file_name) and overwrite != "0":
                x += 1
                file_name = rename_file(directory, csv_file_module_one, "_%i.fasta" % x)

            if len(file_list) > 2500:
                file_name = rename_file(directory, csv_file_module_one, "_temporary.fasta")
                output_name[0] = file_name

            download_fasta(counter_queries, web_env, key, search_index, directory, file_name)

        if accession_taxonomy_output:
            file_name = rename_file(directory, csv_file_module_one, "_taxonomy.tsv")

            if os.path.exists(file_name):
                x = 0
                overwrite = str(input(color.YELLOW +
                                  " Exist already a file with the same name, do you want to overwrite it? " +
                                  "(0 = yes, 1 = no) "
                                  + color.END))
            else:
                overwrite = 1

            while os.path.exists(file_name) and overwrite != "0":
                x += 1
                file_name = rename_file(directory, csv_file_module_one, "_taxonomy_%i.tsv" % x)

            if len(file_list) > 2500:
                file_name = rename_file(directory, csv_file_module_one, "_temporary.tsv")
                output_name[1] = file_name
            download_accession_taxonomy(counter_queries, web_env, key, search_index, directory, file_name)

        if enrich_output:
            download_enrich_file(counter_queries, web_env, key, search_index, directory, file_name=csv_file_module_one)

def download_enrich_file(counter_id, webenv,query_key,query,folder_path, file_name=None):

    print("\nStart downloading enriched output...\n")
    batch_size = 200  # Batch_size value limit our queries to NCBI so we don't get blacklisted
    missing_part = 0
    wm_all=[]
    coordinates = []  # list of coordinates with their respective gene and organism
    countries = []  # list of countries whenever we don't find coordinates from the record
    no_geo = 0  # Counting how many records didn't put the coordinates nor the country value
    genes = []  # List of genes for the menu on the world map

    if file_name is not None:
        title_map = file_name
    else:
        title_map = query

    logging.info(' Downloading info for enrich file...')

    if counter_id > 400000:
        print(webenv[0])
        indice = 0
        div_n = counter_id / 200000
        list_counter = []
        if div_n.is_integer():
            for num in range(0, int(div_n)):
                list_counter.append(200000)
        else:
            for num in range(0, int(div_n)):
                list_counter.append(200000)
            list_counter.append(counter_id-200000*int(div_n))

        #print(list_counter)

        for idx, mini_counter in enumerate(list_counter):
            for start in range(0, mini_counter, batch_size):

                logging.info(' Downloading data for enrich file: %i of %i' % (start + 1, mini_counter))

                data = efetch_call(start, "nuccore", mini_counter,webenv[indice],
                                   query_key, "gbc")

                if data is None:
                    missing_part = 1
                    print("\nMissing part!")
                    continue

                else:

                    try:
                        root = ET.fromstring(data)
                        INSDSeqs = root.findall('.//INSDSeq')

                    except ET.ParseError as PE:
                        print("\nError while parsing at %i" % start)
                        print("Error: %s" % PE)
                        continue

                    for INSDSeq in INSDSeqs:  # Iter every entry of the INSDSeq xml file
                        INSDQualifiers = INSDSeq.findall('.//INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier')
                        lat_lon = None
                        lat = "NA"
                        lon = "NA"
                        country = "NA"
                        gene = []
                        accession = "NA"
                        organism = "NA"

                        for INSDSeq_tags in INSDSeq:
                            if INSDSeq_tags.tag == "INSDSeq_primary-accession":
                                accession = str(INSDSeq_tags.text)

                        for INSDQualifier in INSDQualifiers:

                            if INSDQualifier[0].text == 'organism':
                                organism = str(INSDQualifier[1].text)

                            if INSDQualifier[0].text == 'lat_lon':
                                lat_lon = INSDQualifier[1].text
                                while(True):
                                    try:
                                        if "S" in lat_lon:
                                            lat = 0 - float(lat_lon[0: lat_lon.index("S") - 1])

                                            if "W" in lat_lon:
                                                lon = 0 - float(lat_lon[lat_lon.index("S") + 1: lat_lon.index("W") - 1])
                                                break
                                            else:
                                                lon = float(lat_lon[lat_lon.index("S") + 1: lat_lon.index("E") - 1])
                                                break
                                        else:
                                            lat = float(lat_lon[0: lat_lon.index("N") - 1])

                                            if "W" in lat_lon:
                                                lon = 0 - float(lat_lon[lat_lon.index("N") + 1: lat_lon.index("W") - 1])
                                                break
                                            else:
                                                lon = float(lat_lon[lat_lon.index("N") + 1: lat_lon.index("E") - 1])
                                                break
                                    except ValueError as err:
                                        break

                            if INSDQualifier[0].text == 'country':
                                country = str(INSDQualifier[1].text)
                                if ":" in country:
                                    country = country.split(':')[0]

                            if INSDQualifier[0].text == 'gene':
                                gene = INSDQualifier[1].text
                                genes.append(gene.replace(":", " "))

                        wm_all.append({'accession':accession,'org': organism, 'country': country, 'lat': lat, 'lon': lon, 'gene': gene})

                        if lat_lon is not None:
                            coordinates.append({'org': organism, 'lat': lat, 'lon': lon, 'gene': gene})
                            genes.append(gene)

                        elif country != "NA":
                            countries.append({'org': organism, 'country': country, 'gene': gene})

                        else:
                            no_geo += 1
            indice += 1
    else:
        for start in range(0, counter_id, batch_size):

            logging.info(' Downloading data for enrich file: %i of %i' % (start + 1, counter_id))

            data = efetch_call(start, "nuccore", counter_id, webenv, query_key, "gbc")

            if data is None:
                missing_part = 1
                print("\nMissing part!")
                continue

            else:

                try:
                    root = ET.fromstring(data)
                    INSDSeqs = root.findall('.//INSDSeq')

                except ET.ParseError as PE:
                    print("\nError while parsing at %i" % start)
                    print("Error: %s" % PE)
                    continue

                for INSDSeq in INSDSeqs:  # Iter every entry of the INSDSeq xml file
                    INSDQualifiers = INSDSeq.findall('.//INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier')
                    lat_lon = None
                    lat = "NA"
                    lon = "NA"
                    country = np.nan
                    gene = []
                    accession = "NA"
                    organism = "NA"

                    for INSDSeq_tags in INSDSeq:
                        if INSDSeq_tags.tag == "INSDSeq_primary-accession":
                            accession = str(INSDSeq_tags.text)

                    for INSDQualifier in INSDQualifiers:

                        if INSDQualifier[0].text == 'organism':
                            organism = str(INSDQualifier[1].text)

                        if INSDQualifier[0].text == 'lat_lon':
                            lat_lon = INSDQualifier[1].text
                            while(True):
                                try:
                                    if "S" in lat_lon:
                                        lat = 0 - float(lat_lon[0: lat_lon.index("S") - 1])

                                        if "W" in lat_lon:
                                            lon = 0 - float(lat_lon[lat_lon.index("S") + 1: lat_lon.index("W") - 1])
                                            break
                                        else:
                                            lon = float(lat_lon[lat_lon.index("S") + 1: lat_lon.index("E") - 1])
                                            break
                                    else:
                                        lat = float(lat_lon[0: lat_lon.index("N") - 1])

                                        if "W" in lat_lon:
                                            lon = 0 - float(lat_lon[lat_lon.index("N") + 1: lat_lon.index("W") - 1])
                                            break
                                        else:
                                            lon = float(lat_lon[lat_lon.index("N") + 1: lat_lon.index("E") - 1])
                                            break
                                except ValueError as err:
                                    break

                        if INSDQualifier[0].text == 'country':
                            country = str(INSDQualifier[1].text)
                            if ":" in country:
                                country = country.split(':')[0]

                        if INSDQualifier[0].text == 'gene':
                            gene = INSDQualifier[1].text
                            genes.append(gene.replace(":", " "))

                    wm_all.append({'accession':accession,'org': organism, 'country': country, 'lat': lat, 'lon': lon, 'gene': gene})

                    if lat_lon is not None:
                        coordinates.append({'org': organism, 'lat': lat, 'lon': lon, 'gene': gene})
                        genes.append(gene)

                    elif country != "NA":
                        countries.append({'org': organism, 'country': country, 'gene': gene})

                    else:
                        no_geo += 1

    df_wm=pd.DataFrame(wm_all)
    df_wm.to_csv(os.path.join(folder_path,title_map+'_enriched.tsv'), sep='\t',index=False)
    update_progress(counter_id,counter_id)
    print("\nEnriched dataframe created!\n")

load_configurations("alberto.brusati@gmail.com","0ae434ddfd0897bdefe6398398d80ad12809")
db_creation(text_search='txid8832',accession_taxonomy_output=False, fasta_output=False, marker_output=False, top10_plot=False , enrich_output=True)
#db_creation(file_search='example/A_accession_list_example.tsv',input_file_type='A',accession_taxonomy_output=True, fasta_output=True, marker_output=True,top10_plot=True)
