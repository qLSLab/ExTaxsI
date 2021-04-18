import time

# Biopython
from Bio import Entrez
from Bio.Entrez import Parser as ps

from urllib.error import HTTPError
import http.client as http

# For log files
import logging

import ETI_Lib as ETIL
import ETI_Lib_Interface as ETILI
import ETI_Lib_Retrival as ETILRET


def database_module(plot_or_not, file_pos, file, choice, output_name):
    # if file_pos is more than 0 means the tool is resuming the download from that position in the file
    if file_pos > 0:
        search_list = [file[0]]

        end = min(len(file), file_pos + 2500)
        for pos in range(file_pos, end):
            search_list.append(file[pos])

        if search_list[0] == "accession":
            search_index = [str(search_list[1])]

            for i in range(2, len(search_list)):
                search_index = "{0} or {1}".format(search_index,
                                                   search_list[i])

        else:
            search_index = "(({0}) AND ({1}))".format(search_list[1],
                                                      search_list[0])

            for i in range(2, len(search_list)):
                search_index = "{2} OR (({0}) AND ({1}))".format(
                    search_list[i], search_list[0], search_index)

        # variable to prevent the stop of esearch if you get a momentary internet problem
        server_error = 0

        while server_error < 11:

            try:
                handles = Entrez.esearch(db="nucleotide",
                                         term=search_index,
                                         usehistory="y",
                                         idtype="acc")
                result = Entrez.read(handles, validate=False)
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
                        print(
                            "Maybe NCBI's servers are too busy, if it still doesn't work, try later"
                        )

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
            print("Connection impossible with NCBI, can't request %s" %
                  search_index)
            logging.warning(
                "Connection impossible with NCBI, can't request %s" %
                search_index)
            database_module(None, pos, file, choice, output_name)

        counter_queries = int(
            result["Count"]
        )  # With counter_query we know how many IDs we found
        # results["count"] contains the results counts found with the query
        # This number is needed to define how many ID's we will search and download
        # Default value is 20 so without it NCBI will let us download only a bunch of data instead of the full list

        if counter_queries == 0:
            print("\nNo results found with: \n", search_index, "\n")
            logging.warning("No results found with: %s" % search_index)
            database_module(None, pos, file, choice, output_name)

        # variable to prevent the stop of esearch if you get a momentary internet problem
        server_error = 0

        while server_error < 11:

            try:
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
                        print(
                            "Maybe NCBI's servers are too busy, if it still doesn't work, try later"
                        )

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
            print("Connection impossible with NCBI, can't request %s" %
                  search_index)
            logging.warning(
                "Connection impossible with NCBI, can't request %s" %
                search_index)

        if (choice == 1 or choice == 4) and counter_queries > 0:
            ext = "_{0}_temporary.fasta".format(file_pos)
            file_name = ETIL.rename_file(directory, output_name[0], ext)

            download_fasta(counter_queries, web_env, key, search_index,
                           directory, file_name)

        if (choice == 2 or choice == 4) and counter_queries > 0:
            ext = "_{0}_temporary.tsv".format(file_pos)
            file_name = rename_file(directory, output_name[1], ext)

            download_accession_taxonomy(counter_queries, web_env, key,
                                        search_index, directory, file_name)

        if end >= len(file):
            if output_name[0] is not None:
                print("Merging all temporary files into %s" % output_name[0])

                with open(output_name[0], 'wb') as outfile:
                    for filename in glob.glob(
                            '{0}*_temporary.fasta'.format(directory)):

                        if filename == output_name[0]:
                            # don't want to copy the output into the output
                            continue

                        with open(filename, 'rb') as readfile:
                            outfile.write(readfile.read())

            if output_name[1] is not None:
                print("Merging all temporary files into %s" % output_name[1])

                with open(output_name[1], 'wb') as outfile:

                    for filename in glob.glob(
                            '{0}*_temporary.tsv'.format(directory)):

                        if filename == output_name[1]:
                            # don't want to copy the output into the output
                            continue

                        with open(filename, 'rb') as readfile:
                            outfile.write(readfile.read())

        else:
            database_module(None, pos, file, choice, output_name)

    # Here is how the module starts normally
    else:
        if plot_or_not is None:  # Creating the menu for standalone module (when not called from statistical module)
            print("\n--- DATABASE MODULE --- \n")

        csv_manual = str(
            input(
                "Do you want to use a manual search or through a file? (f > file, m > manual): "
            ))

        while len(csv_manual) <= 0 or csv_manual[0] not in ("f", "F", "m",
                                                            "M"):
            print(
                "Wrong input, write <f> for csv or <m> for manual! (Without <>)\n"
            )
            csv_manual = input(
                "Do you want to use a manual search or through a file? (f > file, m > manual) "
            )

        if csv_manual[0] in (
                "f",
                "F"):  # Considerare la tipologia di file da avere come input

            # Initializing while cycle
            input_list_type = None

            while input_list_type not in ("A", "a", "t", "T", "o", "O"):
                input_list_type = input(
                    "\n'A' if your csv file is a list of Accession\n" +
                    "'T' if it's a list of taxonomy IDs\n" +
                    "'O' if it's a list of organism name\n" +
                    ">> Enter the choice: ")

                if input_list_type in ("T", "t", "O", "o"):
                    # POSSIBLE INPUT -> GENE IN THE COLUMN 0
                    file_list = [
                        str(
                            input(
                                "\nEnter <0> for all genes\nFor multiple genes use 'OR' example:\n- COX1 OR 16S\n"
                                "- COI or COX1 or CO1\n- 0\n>> Enter the gene name:"
                            ))
                    ]

                elif input_list_type in ("A", "a"):
                    file_list = ["accession"]

                else:
                    print("\nEnter a valid choice! \n")

            while True:
                csv_file_module_one = input(
                    "\n>> Enter your csv or tsv file path: ")
                print("")

                try:
                    logging.info("Loading csv file...")
                    with open(csv_file_module_one, 'r') as r:

                        while csv_file_module_one[-3:] not in ('tsv', 'csv'):
                            print(
                                "File extension not CSV or TSV, please use the right format."
                            )
                            csv_file_module_one = input(
                                "Enter your csv or tsv file path: ")

                        if csv_file_module_one[-3:] == 'csv':
                            file_csv_query = csv.reader(r, delimiter=',')

                        if csv_file_module_one[-3:] == 'tsv':
                            file_csv_query = csv.reader(r, delimiter='\t')

                        for line in file_csv_query:
                            file_list.append(line[0])

                        logging.info("Loading completed!")
                        break

                except IOError as err:
                    logging.warning(
                        'No file found or permission denied, please check your file or location: %s'
                        % err)
                    print(
                        "No file found or permission denied, please check your file or location"
                    )
                    print("File location: ", csv_file_module_one, "\n")

            while "/" in csv_file_module_one:
                csv_file_module_one = csv_file_module_one[csv_file_module_one.
                                                          index("/") + 1:]

            if plot_or_not is None:
                print('Choose one of the following options: \n ',
                      '1: Download FASTA files \n ',
                      '2: Download files with accession and taxonomy \n ',
                      '3: Retrieve DNA markers and genes \n ',
                      '4: All the above options \n ',
                      '5: Make a new query \n ', '6: Main menu  \n')

                while True:

                    try:
                        menu_choice = int(input('Enter the NUMBER: '))

                        if menu_choice not in (1, 2, 3, 4, 5, 6):
                            print(
                                "Wrong input, please write only the choice number. \n"
                            )
                            continue

                        break

                    except ValueError:
                        print(
                            "Error, insert only the numeric value of your choice"
                        )
                        continue

                if menu_choice == 5:
                    clear()
                    database_module(None, 0, None, None, [None, None])

                if menu_choice == 6:
                    clear()
                    main_menu()

                if menu_choice == 3 or menu_choice == 4:
                    first_element = True

                    if create_folder(
                            "./download/{0}/".format(csv_file_module_one)):
                        marker_folder = "./download/{0}/".format(
                            csv_file_module_one)

                    else:
                        marker_folder = directory

                    for element in file_list:

                        if first_element:
                            first_element = False
                            continue

                        if input_list_type in ("o",
                                               "O"):  # File list of Organism

                            if file_list[0] == "0":
                                search_index = element

                            else:
                                search_index = "({0} AND {1})".format(
                                    element, file_list[0])

                        elif input_list_type in ("t",
                                                 "T"):  # File list of TaxIDs

                            if file_list[0] == "0":
                                search_index = "{1}{0}{2}".format(
                                    element, "txid", "[Organism:noexp]")

                            else:
                                search_index = "({2} AND {1}{0}{3})".format(
                                    element, "txid", file_list[0],
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
                                    logging.warning(
                                        "No results found with: %s" %
                                        search_index)
                                    break

                                # results["count"] contains the results counts found with the query
                                # This number is needed to define how many ID's we will search and download
                                # Default value is 20 and NCBI will let us download only 20 ids instead of the full list
                                handles = Entrez.esearch(
                                    db="nucleotide",
                                    term=search_index,
                                    retmax=counter_queries)

                                result = Entrez.read(handles, validate=False)
                                id_list = result["IdList"]
                                epost_xml = Entrez.epost("nucleotide",
                                                         id=",".join(id_list))
                                result = Entrez.read(epost_xml, validate=False)
                                web_env = result["WebEnv"]
                                key = result["QueryKey"]
                                break

                            except HTTPError as err:
                                server_error += 1

                                if 500 <= err.code <= 599:
                                    logging.warning(
                                        ' Connection error to the server %s' %
                                        err)
                                    print(
                                        '\nConnection error to the server %s' %
                                        err)
                                    print(search_index)
                                    time.sleep(60)
                                    logging.warning(" Attempt %i of 10" %
                                                    server_error)
                                    print("Attempt %i of 10" % server_error)
                                    time.sleep(10)
                                    if server_error > 6:
                                        print(
                                            "Maybe NCBI's servers are too busy, if it still doesn't work, try later"
                                        )

                                else:
                                    print("\nError in the connection: %s" %
                                          err)
                                    logging.error("Exception occurred",
                                                  exc_info=True)
                                    print("Attempt %i of 10" % server_error)
                                    time.sleep(20)

                            except OSError as err:
                                print("\nError in the connection: %s" % err)
                                logging.error("Exception occurred",
                                              exc_info=True)
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
                            print(
                                "Connection impossible with NCBI, can't request %s"
                                % search_index)
                            logging.warning(
                                "Connection impossible with NCBI, can't request %s"
                                % search_index)
                            #Bloccare?

                        file_name = rename_file(marker_folder, element,
                                                "_gene_list.tsv")
                        print(file_name)

                        if os.path.exists(file_name):
                            x = 0
                            overwrite = str(
                                input(
                                    ETIL.color.YELLOW +
                                    "Exist already a file with the same name, do you want to overwrite it? "
                                    + "(0 = yes, 1 = no) " + ETIL.color.END))

                        else:
                            overwrite = 1

                        while os.path.exists(file_name) and overwrite != "0":
                            x += 1
                            file_name = rename_file(directory, element,
                                                    "_gene_list_%i.tsv" % x)

                        if counter_queries != 0:
                            download_gene_markers(counter_queries, web_env,
                                                  key, search_index, directory,
                                                  file_name)
                            time.sleep(2)

                    merge_gene_top10(None, marker_folder)

                if len(
                        file_list
                ) > 2500:  # NCBI rejects more than 2500 - 2700 terms in a unique call so we will split them
                    search_list = [file_list[0]]
                    for pos in range(1, 2500):
                        search_list.append(file_list[pos])
                else:
                    search_list = file_list

                if input_list_type in ("o", "O"):  # File list of Organism
                    if search_list[0] == "0":
                        search_index = "{0} OR {1}".format(
                            search_list[1], search_list[2])
                        for i in range(3, len(search_list)):
                            search_index = "{0} OR {1}".format(
                                search_index, search_list[i])

                    else:
                        search_index = "({0} AND {1})".format(
                            search_list[1], search_list[0])

                        for i in range(2, len(search_list)):
                            search_index = "{2} OR ({0} AND {1})".format(
                                search_list[i], search_list[0], search_index)

                elif input_list_type in ("t", "T"):  # File list of TaxIDs

                    if search_list[0] == "0":
                        search_index = "{2}{0}{3} OR {2}{1}{3}".format(
                            search_list[1], search_list[2], "txid",
                            "[Organism:noexp]")
                        for i in range(3, len(search_list)):
                            search_index = "{0} OR {2}{1}{3}".format(
                                search_index, search_list[i], "txid",
                                "[Organism:noexp]")

                    else:
                        search_index = "({2} AND {1}{0}{3})".format(
                            search_list[1], "txid", search_list[0],
                            "[Organism:noexp]")
                        for i in range(2, len(search_list)):
                            search_index = "{0} OR ({3} AND {2}{1}{4})".format(
                                search_index, search_list[i], "txid",
                                search_list[0], "[Organism:noexp]")

                else:  # File list of ACCESSIONS
                    search_index = [str(search_list[1])]

                    for i in range(2, len(search_list)):
                        search_index = "{0} or {1}".format(
                            search_index, search_list[i])

                # variable to prevent the stop of esearch if you get a momentary internet problem
                server_error = 0

                while server_error < 11:

                    try:
                        handles = Entrez.esearch(db="nucleotide",
                                                 term=search_index,
                                                 usehistory="y",
                                                 idtype="acc")
                        break

                    except HTTPError as err:
                        server_error += 1

                        if 500 <= err.code <= 599:
                            logging.warning(
                                ' Connection error to the server %s' % err)
                            print('\nConnection error to the server %s' % err)
                            print(search_index)
                            time.sleep(60)
                            logging.warning(" Attempt %i of 10" % server_error)
                            print("Attempt %i of 10" % server_error)
                            time.sleep(10)
                            if server_error > 6:
                                print(
                                    "Maybe NCBI's servers are too busy, if it still doesn't work, try later"
                                )

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

                    except http.client.IncompleteRead as err:
                        print(err)
                        print("\n Retrying in a few seconds")
                        time.sleep(10)

                    except RuntimeError as err:
                        print(err)
                        print("\n Retrying in a few seconds")
                        time.sleep(10)

                if server_error >= 11:
                    print("Connection impossible with NCBI, can't request %s" %
                          search_index)
                    logging.warning(
                        "Connection impossible with NCBI, can't request %s" %
                        search_index)
                    #Bloccare?

                result = Entrez.read(handles, validate=False)
                counter_queries = int(
                    result["Count"]
                )  # With counter_query we know how many IDs we found

                if counter_queries == 0:
                    print("\nNo results found with: \n", search_index, "\n")
                    logging.warning("No results found with: %s" % search_index)

                # results["count"] contains the results counts found with the query
                # This number is needed to define how many ID's we will search and download
                # Default value is 20 so without it NCBI will let us download only a bunch of data instead of the full list
                handles = Entrez.esearch(db="nucleotide",
                                         term=search_index,
                                         retmax=counter_queries)

                result = Entrez.read(handles, validate=False)
                id_list = result["IdList"]
                epost_xml = Entrez.epost("nucleotide", id=",".join(id_list))
                result = Entrez.read(epost_xml, validate=False)
                web_env = result["WebEnv"]
                key = result["QueryKey"]

                if menu_choice == 1 or menu_choice == 4:
                    file_name = rename_file(directory, csv_file_module_one,
                                            ".fasta")

                    if os.path.exists(file_name):
                        x = 0
                        overwrite = str(
                            input(
                                ETIL.color.YELLOW + "0 => yes, 1 => no" +
                                "Exist already a file with the same name, do you want to overwrite it? "
                                + ETIL.color.END))
                    else:
                        overwrite = 1

                    while os.path.exists(file_name) and overwrite != "0":
                        x += 1
                        file_name = rename_file(directory, csv_file_module_one,
                                                "_%i.fasta" % x)

                    if len(file_list) > 2500:
                        file_name = rename_file(directory, csv_file_module_one,
                                                "_temporary.fasta")
                        output_name[0] = file_name

                    download_fasta(counter_queries, web_env, key, search_index,
                                   directory, file_name)

                if menu_choice == 2 or menu_choice == 4:
                    file_name = rename_file(directory, csv_file_module_one,
                                            "_taxonomy.tsv")

                    if os.path.exists(file_name):
                        x = 0
                        overwrite = str(
                            input(
                                ETIL.color.YELLOW +
                                " Exist already a file with the same name, do you want to overwrite it? "
                                + "(0 = yes, 1 = no) " + ETIL.color.END))
                    else:
                        overwrite = 1

                    while os.path.exists(file_name) and overwrite != "0":
                        x += 1
                        file_name = rename_file(directory, csv_file_module_one,
                                                "_taxonomy_%i.tsv" % x)

                    if len(file_list) > 2500:
                        file_name = rename_file(directory, csv_file_module_one,
                                                "_temporary.tsv")
                        output_name[1] = file_name
                    download_accession_taxonomy(counter_queries, web_env, key,
                                                search_index, directory,
                                                file_name)

                if len(file_list) > 2500:
                    database_module(None, 2500, file_list, menu_choice,
                                    output_name)

                if input("Do you want to go to main menu or stop?\n"
                         ">> 0 for main menu, anything else to close: ") == 0:
                    main_menu()

            else:
                if len(
                        file_list
                ) > 2500:  # NCBI rejects more than 2500 - 2700 terms in a unique call so we will split them
                    search_list = [file_list[0]]
                    for pos in range(1, 2500):
                        search_list.append(file_list[pos])
                else:
                    search_list = file_list

                if input_list_type in ("o", "O"):  # File list of Organism
                    if search_list[0] == "0":
                        search_index = "{0} OR {1}".format(
                            search_list[1], search_list[2])
                        for i in range(3, len(search_list)):
                            search_index = "{0} OR {1}".format(
                                search_index, search_list[i])

                    else:
                        search_index = "({0} AND {1})".format(
                            search_list[1], search_list[0])

                        for i in range(2, len(search_list)):
                            search_index = "{2} OR ({0} AND {1})".format(
                                search_list[i], search_list[0], search_index)

                elif input_list_type in ("t", "T"):  # File list of TaxIDs

                    if search_list[0] == "0":
                        search_index = "{2}{0}{3} OR {2}{1}{3}".format(
                            search_list[1], search_list[2], "txid",
                            "[Organism:noexp]")
                        for i in range(3, len(search_list)):
                            search_index = "{0} OR {2}{1}{3}".format(
                                search_index, search_list[i], "txid",
                                "[Organism:noexp]")

                    else:
                        search_index = "({2} AND {1}{0}{3})".format(
                            search_list[1], "txid", search_list[0],
                            "[Organism:noexp]")
                        for i in range(2, len(search_list)):
                            search_index = "{0} OR ({3} AND {2}{1}{4})".format(
                                search_index, search_list[i], "txid",
                                search_list[0], "[Organism:noexp]")

                else:  # File list of ACCESSIONS
                    search_index = [str(search_list[1])]

                    for i in range(2, len(search_list)):
                        search_index = "{0} or {1}".format(
                            search_index, search_list[i])

                # variable to prevent the stop of esearch if you get a momentary internet problem
                server_error = 0

                while server_error < 11:

                    try:
                        handles = Entrez.esearch(db="nucleotide",
                                                 term=search_index,
                                                 usehistory="y",
                                                 idtype="acc")
                        break

                    except HTTPError as err:
                        server_error += 1

                        if 500 <= err.code <= 599:
                            logging.warning(
                                ' Connection error to the server %s' % err)
                            print('\nConnection error to the server %s' % err)
                            print(search_index)
                            time.sleep(60)
                            logging.warning(" Attempt %i of 10" % server_error)
                            print("Attempt %i of 10" % server_error)
                            time.sleep(10)
                            if server_error > 6:
                                print(
                                    "Maybe NCBI's servers are too busy, if it still doesn't work, try later"
                                )

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
                    print("Connection impossible with NCBI, can't request %s" %
                          search_index)
                    logging.warning(
                        "Connection impossible with NCBI, can't request %s" %
                        search_index)
                    #Bloccare?

                result = Entrez.read(handles, validate=False)
                counter_queries = int(
                    result["Count"]
                )  # With counter_query we know how many IDs we found

                if counter_queries == 0:
                    print("\nNo results found with: \n", search_index, "\n")
                    logging.warning("No results found with: %s" % search_index)

                # results["count"] contains the results counts found with the query
                # This number is needed to define how many ID's we will search and download
                # Default value is 20 so without it NCBI will let us download only a bunch of data instead of the full list
                handles = Entrez.esearch(db="nucleotide",
                                         term=search_index,
                                         retmax=counter_queries)

                result = Entrez.read(handles, validate=False)
                id_list = result["IdList"]
                epost_xml = Entrez.epost("nucleotide", id=",".join(id_list))
                result = Entrez.read(epost_xml, validate=False)
                web_env = result["WebEnv"]
                key = result["QueryKey"]

                if plot_or_not == "scatter":
                    return download_accession_taxonomy(counter_queries,
                                                       web_env, key,
                                                       search_index, directory,
                                                       None)

                else:
                    return {
                        'counter_id': counter_queries,
                        'webenv': web_env,
                        'query_key': key,
                        'query': search_index,
                        'folder_path': directory,
                        'file_input': None
                    }

    # Otherwise search_term is a string (manual mode)
        else:
            print(
                "|----------------------------------------------------------------------------|"
            )
            print(
                "|---------------------------  SEARCH FIELD TAGS  ----------------------------|"
            )
            print(
                "|   Accession = [ACCN]     All Fields      = [ALL]      Gene Name = [GENE]   |"
            )
            print(
                "|   Organism  = [ORGN]     Sequence Length = [SLEN]     Filter    = [FILT]   |"
            )
            print(
                "|----------------------------------------------------------------------------|\n"
            )
            print(
                "Examples of queries:\n"
                "1) odonata[ORGN] AND (latest[filter] AND all[FILT] NOT anomalous[filter]) "
            )
            print(
                "2) txid33208[Organism:exp] AND (18S OR SSU) NOT (mitochondrial OR "
                + "complete genome OR whole genome)\n")

            # Possible HELP function?
            counter_query = 0

            while counter_query == 0:
                server_errors = 0

                search_term = input(">> Enter your query: ")

                while server_errors < 11:
                    server_errors += 1

                    try:
                        handle = Entrez.esearch(db="nucleotide",
                                                term=search_term,
                                                usehistory="y",
                                                idtype="acc")
                        break

                    except HTTPError as err:
                        server_errors += 1

                        if 500 <= err.code <= 599:
                            logging.warning(
                                ' Connection error to the server %s' % err)
                            print('Connection error to the server %s' % err)
                            logging.warning(" Attempt %i of 10" %
                                            server_errors)
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
                    print(
                        "Possible maintenance at NCBI servers if the problem is not in your connection, "
                        + "try again some minutes later")
                    logging.warning(
                        "Connection impossible with NCBI, can't request %s" %
                        search_term)
                    continue

                results = Entrez.read(handle, validate=False)

                # results["count"] contains the results counts found with the query
                # This number is needed to define how many ID's we will search and download
                # Default value is 20 means NCBI will let us download only a bunch of data instead of the full list
                counter_query = int(results["Count"])

                if counter_query == 0:
                    print(
                        "\nNo results found with this query: %s\nMaybe a typo?\n"
                        % search_term)
                    continue
                else:
                    print(
                        "\nThe number of sequence found is: {0}{2}{1}".format(
                            ETILI.color.GREEN, ETILI.color.END, counter_query))
                    if counter_query > 1000000:
                        print(ETILI.color.RED + "High number of records!" +
                              ETILI.color.END +
                              " It can take a while to load..")

                server_errors = 0

                while server_errors < 8:
                    server_errors += 1

                    try:
                        handle = Entrez.esearch(db="nucleotide",
                                                term=search_term,
                                                retmax=counter_query)
                        results = Entrez.read(handle, validate=False)

                        id_list = results["IdList"]

                        if len(id_list) > 1000000:
                            webenv = []
                            query_key = []

                            for index in range(0, len(id_list), 1000000):
                                ###
                                print("index")
                                print(index)
                                ###
                                mini_id_list = []
                                end = min(index + 1000000, len(id_list))

                                for id in range(index, end):
                                    mini_id_list.append(id_list[id])

                                # Using epost we will lighten the burden on Entrez
                                post_xml = Entrez.epost(
                                    "nucleotide", id=",".join(mini_id_list))
                                results = Entrez.read(post_xml, validate=False)

                                # Variable needed to use ncbi history feature
                                webenv.append(results["WebEnv"])
                                query_key.append(results["QueryKey"])

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
                            logging.warning(
                                ' Connection error to the server %s' % err)
                            print('\nConnection error to the server %s' % err)
                            logging.warning(" Attempt %i of 7" % server_errors)
                            print("Attempt %i of 7" % server_errors)
                            time.sleep(10)
                            if server_errors > 6:
                                print(
                                    "Maybe NCBI's servers are too busy, if it still doesn't work, try later"
                                )

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
                    print("Connection impossible with NCBI, can't request %s" %
                          search_term)
                    logging.warning(
                        "Connection impossible with NCBI, can't request %s" %
                        search_term)
                    continue

            if plot_or_not is None:
                while True:
                    print('Choose one of the following options: \n ',
                          '1: Download FASTA files \n ',
                          '2: Download files with accession and taxonomy \n ',
                          '3: Retrieve DNA markers and genes \n ',
                          '4: All the above options \n ',
                          '5: Make a new query \n ', '6: Main menu \n ')

                    choice = None
                    while choice is None:
                        try:
                            choice = int(input('Enter the NUMBER: '))

                            if choice not in (1, 2, 3, 4, 5, 6):
                                print(
                                    "Wrong input, please enter the number of your choice. \n"
                                )
                                choice = None
                                continue

                        except ValueError:
                            print(
                                "Error, insert only the numeric value of your choice"
                            )
                            continue

                    if choice == 1 or choice == 4:
                        ETILRET.download_fasta(counter_query, webenv,
                                               query_key, search_term,
                                               directory, None)

                    if choice == 2 or choice == 4:
                        ETILRET.download_accession_taxonomy(
                            counter_query, webenv, query_key, search_term,
                            directory, None)

                    if choice == 3 or choice == 4:
                        ETILRET.download_gene_markers(
                            counter_query,
                            webenv,
                            query_key,
                            search_term,
                            directory,
                            None,
                        )

                        ETILRET.merge_gene_top10(search_term, directory)

                    if choice == 5:
                        ETILI.clear()
                        database_module(None, 0, None, None, [None, None])

                    if choice == 6:
                        ETILI.clear()
                        ETILI.main_menu()
                        return

            elif plot_or_not == "scatter":
                return ETILRET.download_accession_taxonomy(
                    counter_query, webenv, query_key, search_term, directory,
                    None)

            elif plot_or_not == "world":
                return {
                    'counter_id': counter_query,
                    'webenv': webenv,
                    'query_key': query_key,
                    'query': search_term,
                    'folder_path': directory,
                    'file_input': None
                }

        print("See you later!")
