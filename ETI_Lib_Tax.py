import logging
import csv

import ETI_Lib_Interface as ETILI


def taxonomyID_module():
    print("---- TAXONOMY ID CONVERTER MODULE ----\n")
    id_or_name = 0
    file_manual = str(
        input(
            "Do you want to convert taxonomy IDs through a file or manual input? (f > file, m > manual) "
        ))

    while len(file_manual) <= 0 or file_manual[0] not in ("f", "F", "m", "M"):
        print(
            "Wrong input, write <f> for csv or <m> for manual! (Without <>)\n")
        file_manual = input(
            "Do you want to convert taxonomy IDs through a tsv or manual input? (f > file, m > manual) "
        )

    if file_manual[0] in (
            "m", "M"):  # Considerare la tipologia di file da avere come input
        while (True):
            try:
                id_or_name = int(
                    input('1 - To convert TaxIDs to Taxonomy names\n' +
                          '2 - To convert Taxonomy names to TaxIDs\n' +
                          '>> Enter the option number: '))
                if id_or_name == 1:
                    taxa_input = str(
                        input(
                            "\nIf you want to search for multiple taxonomy IDs add '+' as separator \n"
                            + "Example: 6224+458+26351+6223\n" +
                            "\n>> Enter one or more taxonomy IDs: "))
                    taxa_ids = taxa_input.split("+")
                    skip_first_row = False
                    break

                if id_or_name == 2:
                    taxa_input = str(
                        input(
                            "\nIf you want to search for multiple taxonomies name add '+' as separator \n"
                            +
                            "Example: odonata+bufo bufo+arachnida+Hymenoptera\n"
                            + "\n>> Enter one or more taxonomies name: "))
                    taxa_ids = taxa_input.split("+")
                    skip_first_row = False
                    break

            except ValueError as err:
                print('\nPlease retry! Enter 1 or 2!')
                continue

    else:

        while True:
            taxa_id_path = input(">> Enter your taxonomy ID file path: ")
            print("")
            taxa_ids = []

            try:
                logging.info("Loading tsv file...")
                with open(taxa_id_path, 'r') as r:
                    file_csv_query = csv.reader(r,
                                                skipinitialspace=True,
                                                delimiter='\t')

                    for line in file_csv_query:
                        taxa_ids.append(line[0])

                    logging.info("Loading completed!")
                    break

            except IOError as err:
                logging.warning(
                    'No file found or permission denied, please check your file or location: %s'
                    % err)
                print(
                    "No file found or permission denied, please check your file or location"
                )
                print("File location: ", taxa_id_path, "\n")

        file_name = str(input(">> Enter your output file name: "))

        if len(file_name) <= 0:
            print(
                "No name submitted or unexpected symbols, using your input file name + '_taxonomy_ID_output.txt')"
            )
            file_name = "{0}_taxonomy_ID_output.txt".format(taxa_id_path)

        else:
            file_name = "./download/{0}".format(file_name)

        out_handle = open(file_name, "w")

        ## POSSIBILE MENU PER SCEGLIERE CHE LIVELLO TASSONOMICO POSSA SERVIRE
        print(taxa_ids[0])
        if input("Is this row a title? (Y > yes, N > no) ") in ("Y", "y"):
            skip_first_row = True
            print("Ok, i'll skip it")

        else:
            print("Good!\n")
            skip_first_row = False

    if id_or_name == 2:
        print(ncbi.get_name_translator(taxa_ids))

    else:
        for counter, taxa in enumerate(taxa_ids):

            if file_manual[0] in ("m", "M"):
                print("Tax ID: %s" % taxa)
            if file_manual[0] not in ("m", "M"):
                ETILI.update_progress(counter, len(taxa_ids))

            if skip_first_row:
                skip_first_row = False
                continue

            try:
                # Getting entire taxonomy from its taxonomy id
                lineage = ncbi.get_lineage(int(taxa))

            except ValueError:  # In case its id is not found
                logging.error(
                    " There is no taxonomy for %s, ncbi.get_lineage exception"
                    % taxa,
                    exc_info=True)
                if file_manual[0] in ("m", "M"):
                    print(
                        "There is no taxonomy for %s, ncbi.get_lineage exception\n"
                        % taxa)

                else:
                    taxonomy_id_out = "{0}\t{1}{2}".format(
                        taxa, "NA;NA;NA;NA;NA;NA", "\n")
                    out_handle.write(taxonomy_id_out)
                continue

            phylum = clas = order = family = genus = specie = "NA"  # Initializing

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

                    if "species" == lineage_rank[lineage[z]]:
                        rank_tmp = ncbi.get_taxid_translator([lineage[z]])
                        specie = rank_tmp[lineage[z]]

                taxa_list = ";".join(
                    [phylum, clas, order, family, genus, specie])
                if file_manual[0] in ("m", "M"):
                    print("Phylum;Class;Order;Family;Genus;Species\n%s\n" %
                          taxa_list)

                else:
                    taxa_list = ";".join(
                        [phylum, clas, order, family, genus, specie])
                    taxonomy_id_out = "{0}\t{1}{2}".format(
                        taxa, taxa_list, "\n")
                    out_handle.write(taxonomy_id_out)

            else:
                if file_manual[0] in ("m", "M"):
                    print("Not found in Lineage %i" % taxa)

                else:
                    taxa_list = ";".join(
                        [phylum, clas, order, family, genus, specie])
                    taxonomy_id_out = "{0}\t{1}{2}".format(
                        taxa, taxa_list, "\n")
                    out_handle.write(taxonomy_id_out)

        if file_manual[0] not in ("m", "M"):
            ETILI.update_progress(len(taxa_ids), len(taxa_ids))
            out_handle.close()

    print('\nChoose one of the following options: \n ',
          '1 - Repeat with another file or manual input \n ',
          '2 - Main menu \n ', 'CTRL + C: Close ExTaxsI \n ')

    while True:

        try:
            menu_choice = int(input('>> Enter the NUMBER: '))

            if menu_choice not in (1, 2, 3):
                print("Wrong input, please write only the choice number. \n")
                continue

            break

        except ValueError:
            print("Error, insert only the numeric value of your choice")
            continue

    if menu_choice == 1:
        ETILI.clear()
        taxonomyID_module()

    if menu_choice == 2:
        ETILI.clear()
        ETILI.main_menu()

