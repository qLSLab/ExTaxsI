import time

# For log files
import logging

#
import http.client as http
from urllib.error import HTTPError

# Biopython
from Bio import Entrez
from Bio.Entrez import Parser as ps

import ETI_Lib as ETIL
import ETI_Lib_Interface as ETILI


def efetch_call(start, db, counter_id, webenv, query_key, rettype):

    batch_size = 200
    ETILI.update_progress(start, counter_id)  # Download animation bar function

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
                print(ETILI.color.RED + '\nConnection error: %s' % err +
                      ETILI.color.END)
                logging.warning(" Attempt %i of 10" % server_errors)
                print("Attempt %i of 10" % server_errors)
                time.sleep(10)
                if server_errors > 6:
                    print(
                        "Maybe NCBI's servers are too busy, if it still doesn't work, try later"
                    )

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
def download_fasta(counter_id, webenv, query_key, query, folder_path,
                   file_input):
    batch_size = 200  # Batch_size value limit our queries to NCBI so we don't get blacklisted
    if file_input is None:
        file_name = ETIL.rename_file(folder_path, query,
                                     ".fasta")  # Formatting the file name
        out_handle = open(
            file_name,
            "w")  # Creating the fasta file, where we will save our downloads

    logging.info('Downloading Fasta files')
    missing_part = 0

    print("Downloading fasta sequence")
    for start in range(0, counter_id, batch_size):
        logging.info('Downloading fasta sequence: %i of %i' %
                     (start + 1, counter_id))
        data = efetch_call(start, "nucleotide", counter_id, webenv, query_key,
                           "fasta")

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

    ETILI.update_progress(counter_id, counter_id)

    if missing_part == 0:
        logging.info(' Fasta file has been created')
        print("\n ---- Fasta file has been created. ----\n\n")

    else:
        print(
            "\n ---- Fasta file has been created but some errors occurred. ----"
        )
        print(" ---- Check log file for more details. ----\n\n")


# Create a file with ACCESSION iD and six level of Taxonomy (Phylum,Class,Order,Family,Genus,Species)
def download_accession_taxonomy(counter_id, webenv, query_key, query,
                                folder_path, file_input):
    logging.info(' Creating accession file')
    if file_input is None:
        file_name = input("Enter file name: ")
        file_name = ETIL.rename_file(folder_path, file_name, "_taxonomy.tsv")
        out_handle = open(file_name, "w")

    batch_size = 200
    missing_part = 0

    print("Downloading accession with taxonomy file")
    for start in range(0, counter_id, batch_size):
        logging.info(' Downloading accession with taxonomy file: %i of %i' %
                     (start + 1, counter_id))

        data = efetch_call(start, "nucleotide", counter_id, webenv, query_key,
                           "gpc")

        if data is None:
            missing_part = 1
            continue

        else:

            for n in range(len(data)):  # For every downloaded entry
                tax_id = None

                # if taxonomy id is not found, it will give NA for every taxonomy level
                try:
                    tax_id = ncbi.get_name_translator(
                        [data[n]["INSDSeq_organism"]])

                    if tax_id:
                        tax_id = tax_id[data[n]["INSDSeq_organism"]][0]

                    else:
                        tax_id = None

                except ValueError:
                    logging.warning(" There is no taxonomy for %s \n" +
                                    " searching through taxID" %
                                    data[n]["INSDSeq_organism"])

                except sqlite3.OperationalError:
                    logging.warning(" There is no taxonomy for %s \n " +
                                    "searching through taxID" %
                                    data[n]["INSDSeq_organism"])

                if tax_id is None:
                    for z in range(
                            len(data[n]["INSDSeq_feature-table"][0]
                                ["INSDFeature_quals"])):
                        # Searching for taxonomy id in the xml type data
                        try:
                            if "taxon" in data[n]["INSDSeq_feature-table"][0][
                                    "INSDFeature_quals"][z][
                                        "INSDQualifier_value"]:
                                tax_id = data[n]["INSDSeq_feature-table"][0][
                                    "INSDFeature_quals"][z][
                                        "INSDQualifier_value"].replace(
                                            "taxon:", "")
                                break

                        except KeyError:
                            continue

                if tax_id is None:
                    logging.error(
                        " There is no taxonomy for %s, ncbi.get_lineage exception"
                        % data[n]["INSDSeq_organism"],
                        exc_info=True)
                    taxa_list = "NA;NA;NA;NA;NA"
                    accession_list = "{0}\t{1};{2}{3}".format(
                        data[n]["INSDSeq_accession-version"], taxa_list,
                        data[n]["INSDSeq_organism"], "\n")
                    if file_input is None:
                        out_handle.write(accession_list)

                    else:
                        with open(file_input, 'a+') as a_writer:
                            a_writer.write(accession_list)
                    continue

                try:
                    lineage = ncbi.get_lineage(
                        tax_id)  # Getting entire taxonomy from its taxonomy id

                except ValueError:  # In case its id is not found
                    logging.error(
                        " There is no taxonomy for %s, ncbi.get_lineage exception"
                        % data[n]["INSDSeq_organism"],
                        exc_info=True)
                    taxa_list = "NA;NA;NA;NA;NA"
                    accession_list = "{0}\t{1};{2}{3}".format(
                        data[n]["INSDSeq_accession-version"], taxa_list,
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
                accession_list = "{0}\t{1};{2}{3}".format(
                    data[n]["INSDSeq_accession-version"], taxa_list,
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

    ETILI.update_progress(counter_id, counter_id)
    logging.info(' Accession with taxonomy file has been created')

    if missing_part == 0:
        print("\n ---- Accession with taxonomy file has been created. ----\n")

    else:
        print(
            "\n ---- Accession with taxonomy file has been created but some errors occurred. ----"
        )
        print(" ---- Check log file for more details. ---- \n")
    if file_input is None:
        return file_name
    else:
        return


def download_gene_markers(counter, webenv, query_key, search_term, path, file):

    if counter == 0:
        return

    batch_size = 200
    genes_raw = []

    if isinstance(webenv, str):
        for start in range(0, counter, batch_size):

            logging.info(' Downloading data for world map: %i of %i' %
                         (start + 1, counter))
            data = efetch_call(start, "nuccore", counter, webenv, query_key,
                               "gbc")

            if data is None:
                print("\nMissing part!")
                continue

            else:
                try:
                    root = ET.fromstring(data)
                    INSDSeqs = root.findall('.//INSDSeq')

                    for INSDSeq in INSDSeqs:  # Iter every entry of our download
                        INSDQualifiers = INSDSeq.findall(
                            './/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier'
                        )

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
                                    no_double_gene.append(
                                        INSDQualifier[1].text)
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
    else:
        for i in range(len(webenv)):

            if i == (len(webenv) - 1):
                end = counter - (i * 1000000)
            else:
                end = 1000000

            for start in range(0, end, batch_size):
                logging.info(' Downloading data for world map: %i of %i' %
                             (start + 1, end))
                data = efetch_call(start, "nuccore", end, webenv[i],
                                   query_key[i], "gbc")

                if data is None:
                    print("\nMissing part!")
                    continue

                else:

                    try:
                        root = ET.fromstring(data)
                        INSDSeqs = root.findall('.//INSDSeq')

                        for INSDSeq in INSDSeqs:  # Iter every entry of our download
                            INSDQualifiers = INSDSeq.findall(
                                './/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier'
                            )
                            no_double_gene = []

                            for INSDQualifier in INSDQualifiers:

                                if INSDQualifier[0].text == 'gene':
                                    if INSDQualifier[1].text in no_double_gene:
                                        continue

                                    else:
                                        no_double_gene.append(
                                            INSDQualifier[1].text)
                                        genes_raw.append(INSDQualifier[1].text)

                    except ET.ParseError as PE:
                        print("Error while parsing at %i \n" % start)
                        print("Error: %s" % PE)
                        root = None
                        INSDSeqs = None
                        print("Error: %s" % PE)
                        root = None
                        INSDSeqs = None
                        INSDQualifiers = None
                        continue

    genes = []

    #  https://www.genecards.org/cgi-bin/carddisp.pl?gene=MT-CYB
    CYTB = [
        "Mitochondrially Encoded Cytochrome B", "Cytochrome B", "MTCYB",
        "CYTB", "Ubiquinol-Cytochrome-C Reductase Complex Cytochrome B",
        "Cytochrome B-C1", "COB"
    ]
    COI = [
        "COXI", "Cytochrome C Oxidase I", "Cytochrome C Oxidase Subunit I",
        "Cytochrome C Oxidase Subunit 1", "COI", "COX1", "MTCO1",
        "Mitochondrially Encoded Cytochrome C Oxidase I",
        "Cytochrome C Oxidase Subunit I", "CO1"
    ]
    COII = [
        "COXII", "Cytochrome C Oxidase II", "Cytochrome C Oxidase Subunit II",
        "Cytochrome C Oxidase Subunit 2", "COX2", "MTCO2",
        "Mitochondrially Encoded Cytochrome C Oxidase II",
        "Cytochrome C Oxidase Subunit II", "COII", "CO2"
    ]
    COIII = [
        "COXIII", "Cytochrome C Oxidase III",
        "Cytochrome C Oxidase Subunit III", "COX3", "MTCO3", "COIII",
        "Cytochrome C Oxidase Subunit 3",
        "Mitochondrially Encoded Cytochrome C Oxidase III", "CO3",
        "Cytochrome C Oxidase Subunit III"
    ]
    ND1 = [
        "MTND1", "NADH1", "NADH Dehydrogenase 1", "Complex I ND1 Subunit",
        "NADH Dehydrogenase Subunit 1", "ND1",
        "NADH Ubiquinone Oxidoreductase Chain 1",
        "Mitochondrially Encoded NADH Dehydrogenase 1", "NAD1"
    ]
    ND2 = [
        "MTND2", "NADH2", "NADH Dehydrogenase 2", "Complex I ND2 Subunit",
        "NADH Dehydrogenase Subunit 2", "ND2",
        "NADH Ubiquinone Oxidoreductase Chain 2",
        "Mitochondrially Encoded NADH Dehydrogenase 2", "NAD2"
    ]
    ND3 = [
        "MTND3", "NADH3", "NADH Dehydrogenase 3", "Complex I ND3 Subunit",
        "NADH Dehydrogenase Subunit 3", "ND3",
        "NADH Ubiquinone Oxidoreductase Chain 3",
        "Mitochondrially Encoded NADH Dehydrogenase 3", "NAD3"
    ]
    ND4 = [
        "MTND4", "NADH4", "NADH Dehydrogenase 4", "Complex I ND4 Subunit",
        "NADH Dehydrogenase Subunit 4", "ND4",
        "NADH Ubiquinone Oxidoreductase Chain 4",
        "Mitochondrially Encoded NADH Dehydrogenase 4", "NAD4"
    ]
    ND5 = [
        "MTND5", "NADH5", "NADH Dehydrogenase 5", "Complex I ND5 Subunit",
        "NADH Dehydrogenase Subunit 5", "ND5",
        "NADH Ubiquinone Oxidoreductase Chain 5",
        "Mitochondrially Encoded NADH Dehydrogenase 5", "NAD5"
    ]
    ND6 = [
        "MTND6", "NADH6", "NADH Dehydrogenase 6", "Complex I ND6 Subunit",
        "NADH Dehydrogenase Subunit 6", "ND6",
        "NADH Ubiquinone Oxidoreductase Chain 6",
        "Mitochondrially Encoded NADH Dehydrogenase 6", "NAD6"
    ]
    ATP6 = [
        "ATP6", "MTATP6",
        "Mitochondrially Encoded ATP Synthase Membrane Subunit 6",
        "ATP Synthase 6", "ATPase6", "ATP Synthase Subunit A",
        "Mitochondrially Encoded ATP Synthase 6", "ATP Synthase F0 Subunit 6"
    ]
    ATP8 = [
        "ATP8", "MTATP8",
        "Mitochondrially Encoded ATP Synthase Membrane Subunit A6L",
        "ATP Synthase 8", "ATPase8", "ATP Synthase Subunit A",
        "Mitochondrially Encoded ATP Synthase 8", "ATP Synthase F0 Subunit 8",
        "A6L"
    ]

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
    ETILI.update_progress(counter, counter)

    if file is None:
        file = "{0}{1}_gene_list.tsv".format(
            path,
            search_term.replace("(", "").replace(")", ""))

    out_handle = open(file, "w")
    out_handle.write("n_records_found\t{0}\n".format(counter))

    for gene in OrderedDict(gene_counter_list.most_common()):
        out_handle.write("{0}\t{1}\n".format(gene, gene_counter_list[gene]))

    out_handle.close()
    print("File created at %s\n" % file)


def merge_gene_top10(search_term, path):
    ###this part will create a unique dataframe from the multiple files retrieved by list-file search:
    while (True):
        try:
            merge_query = input(
                str("\n\nWould you like to merge all retrieved files genes tsv and plot the top-10 (use only if your query was a list-file)?[y/n]"
                    ))
            yesChoice = ['y', 'Y']
            noChoice = ['n', 'N']
            if merge_query[0] in yesChoice:
                final_df = pd.DataFrame()
                os.chdir(path)
                for file in glob.glob(str('*_gene_list.tsv')):
                    df = pd.read_csv(file, sep='\t',
                                     header=None)  # import as dataFrame
                    df_t = df.T  # calculate transpose
                    df_t.columns = df_t.iloc[0]
                    df_t = df_t.drop([0])
                    taxid = file.split('_')[0]  # extract taxid
                    df_t.insert(0, 'Taxid',
                                taxid)  # add taxid column at 0 position
                    final_df = final_df.append(
                        df_t, ignore_index=True,
                        sort=False)  # add row in final_df
                    final_df = final_df.fillna(0)
                final_df.to_csv('final_gene_df.tsv', index=False, sep='\t')

                #final df to plot:
                final_df = final_df.iloc[:, 2:]
                df_t = final_df.T
                df_t['count'] = df_t.sum(axis=1)
                if len(df_t) >= 10:
                    df_t = df_t.head(10).sort_values(by=['count'],
                                                     ascending=True)
                    df_t.reset_index(level=0, inplace=True)
                    df_top10 = df_t[['index', 'count']]
                    graph = top10_graph(' list ', df_top10['count'],
                                        df_top10['index'])
                    graph.write_image('top10_genes.png')
                else:
                    df_t = df_t.sort_values(by=['count'], ascending=True)
                    df_t.reset_index(level=0, inplace=True)
                    df_top = df_t[['index', 'count']]
                    graph = top10_graph(' list ', df_top['count'],
                                        df_top['index'])
                    graph.write_image('top_used_genes.png')
                print('\nfinal_gene_df.tsv and top-plot created\n\n')
                break
            elif merge_query[0] in noChoice:
                plot_query = input(
                    str("\n\nWould you like to plot the top-10 genes?[y/n]"))
                if plot_query[0] in yesChoice:
                    os.chdir(path)
                    dfplot = pd.read_csv(
                        str(search_term.replace("(", "").replace(")", "")) +
                        '_gene_list.tsv',
                        sep='\t',
                        names=['gene', 'count'],
                        skiprows=1)
                    if len(dfplot) >= 10:
                        df_top10 = dfplot.head(10).sort_values(by=['count'])
                        graph = top10_graph(search_term, df_top10['count'],
                                            df_top10['gene'])
                        graph.write_image('top10_genes_' + str(
                            search_term.replace("(", "").replace(")", "")) +
                                          '.png')
                    else:
                        df_top10 = dfplot.sort_values(by=['count'])
                        graph = top10_graph(search_term, df_top10['count'],
                                            df_top10['gene'])
                        graph.write_image('top_used_genes_' + str(
                            search_term.replace("(", "").replace(")", "")) +
                                          '.png')
                    print('\n\nTOP-10 plotted!!!')
                    break
                elif plot_query[0] in noChoice:
                    break
        except IndexError as err:
            print("\nPlease enter y or n.... try again")

