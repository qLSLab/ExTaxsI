import os
import time

import collections
from collections import OrderedDict
from collections import Counter

# For log files
import logging

# Needed to elaborate csv files
import csv
import numpy as np
import pandas as pd

# Needed to search through xml files
import xml.etree.cElementTree as ET

# import 1plotly.plotly as py
import plotly.graph_objects as go
from plotly.offline import plot

import ETI_Lib as ETIL
import ETI_Lib_Interface as ETILI
import ETI_Lib_DB as ETIDB
import ETI_Lib_Retrival as ETILRET


def statistical_module():
    print("--- STATISTICAL MODULE ---\n")

    print('Choose one of the following options: \n ', '1: Scatter plot \n ',
          '2: World map plot \n ', '3: Sunburst-pie plot \n ',
          '4: Back to main menu \n ', 'CTRL + C: Close ExTaxsI \n ')

    plot_menu_choice = 0

    while plot_menu_choice not in (1, 2, 3, 4, 5):
        try:
            plot_menu_choice = int(input(">> Enter your choice: "))

            if plot_menu_choice not in (1, 2, 3, 4, 5):
                print("Wrong choice, try again")

        except ValueError:
            print("Error, insert only the numeric value of your choice")

    if plot_menu_choice == 1:
        ETILI.clear()
        scatterplot()

    elif plot_menu_choice == 2:
        ETILI.clear()
        worldmap_plot()

    elif plot_menu_choice == 3:
        ETILI.clear()
        sunburst_plot()

    elif plot_menu_choice == 4:
        ETILI.clear()
        ETILI.main_menu()


def scatterplot():
    print("Choose which type of input you will use\n"
          "Enter one of the following options number\n"
          "  1 - Taxonomy file from this tool\n"
          "  2 - Doing a query on NCBI and directly plot it\n"
          "  3 - File input from xxx ecc.. for particular settings")

    while True:
        try:
            query_or_path = int(
                input("\n>> Enter your choice (only the number): "))

            while query_or_path not in (1, 2, 3):
                print("Enter a valid choice, try again. \n")
                query_or_path = input(
                    ">> Enter your choice (only the number): ")

            break

        except ValueError:
            print("Please enter only the numeric value")

    # Until they put an existing file path, it will continue to ask
    while True:

        if query_or_path == 2:
            # To download taxonomy file and getting the path
            plot_file_path = ETIDB.database_module("scatter", 0, None, None,
                                                   [None, None])

        else:
            logging.info(" Getting the taxonomy file path ")
            plot_file_path = str(
                input(
                    ">> Enter the taxonomy file path name (only TSV or CSV files): "
                ))

        try:
            if plot_file_path[-3:] == 'csv':
                with open(plot_file_path, 'r') as axe:
                    r = csv.reader(axe, delimiter=",")
                    r_list = list(r)
                    break

            elif plot_file_path[-3:] == 'tsv':
                with open(plot_file_path, 'r') as axe:
                    r = csv.reader(axe, delimiter="\t")
                    r_list = list(r)
                    break
            else:
                print(
                    "File extension not CSV or TSV, please use the right format."
                )

        except IOError:
            logging.warning(
                " No file found or permission denied, path entered: %s" %
                plot_file_path)
            print(
                "No file found or permission denied, please check your file or location"
            )
            print("File location: ", plot_file_path, "\n")

    taxonomy_list = []
    for line in r_list:
        taxonomy_list.append(line[1].split(";"))

    # Sorting the taxonomy file
    for sorter in range(5, 0, -1):
        taxonomy_list.sort(key=lambda li: li[sorter])

    taxonomy_array = np.array(taxonomy_list)

    print(
        "Now enter the minimum results value needed to be plotted for each taxonomy"
    )
    filter_value = None
    while filter_value is None:

        try:
            filter_value = int(input('>> Enter filter value: '))
            if filter_value < 0:
                print("Enter a valid positive number!")
                filter_value = None

        except ValueError:
            print("Enter a number!")

    title_graph = input(">> Enter the title for the plot: ")

    while "/" in plot_file_path:
        x = plot_file_path.find("/")
        plot_file_path = plot_file_path[x + 1:]

    if "." in plot_file_path:
        x = plot_file_path.find(".")
        plot_file_path = plot_file_path[:x]

    plot_file_name = ETILI.rename_file(directory, plot_file_path,
                                       '_scatterplot.html')
    filter_title = ETILI.rename_file(
        directory,
        "{0}_with_less_than_{1}_results".format(plot_file_path,
                                                str(filter_value)), ".txt")
    if os.path.exists(plot_file_name):
        x = 0
        try:
            overwrite = int(
                input(
                    ETILI.color.YELLOW +
                    "Exist already a file with the same name, do you want to overwrite it? "
                    + "(0 = yes, 1 = no)" + ETILI.color.END))

        except ValueError:
            overwrite = 1

    else:
        overwrite = 1

    while os.path.exists(plot_file_name) and overwrite != 0:
        x += 1
        plot_file_name = ETIL.rename_file(directory, plot_file_path,
                                          '_scatterplot_%i.html' % x)
        filter_title = ETIL.rename_file(
            directory,
            "{0}_with_less_than_{1}_results".format(plot_file_path,
                                                    str(filter_value)),
            "_%i.txt" % x)

    # dict of rgb data to color our points based on their parent level taxonomy
    color_list = {}

    # Initializing the list of traces for scatter plot
    trace = []

    taxa_level = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    # Initializing the button list
    button_counter = []

    for col in range(len(taxonomy_list[0])):
        pos_x = 0
        print("\nElaborating {0}{1}{2}".format(ETILI.color.BOLD,
                                               taxa_level[col],
                                               ETILI.color.END))
        logging.info("Elaborating data of column: %i" % int(col))
        taxonomy_counter = collections.Counter(taxonomy_array[:, col])
        button_counter.append(len(taxonomy_counter.keys()))

        if col > 0:

            if "NA" in taxonomy_array[:, col]:
                # Searching for every NA in list and creating a list of their index
                na_index = [
                    i for i, v in enumerate(taxonomy_array[:, col])
                    if v == 'NA'
                ]
                # Unique list of their parent taxa with NA as son
                unique_parent_taxa_na = list(
                    np.unique(taxonomy_array[na_index, col - 1], axis=0))

            # Unique parent taxonomy
            unique_parent_taxa = list(
                np.unique(taxonomy_array[:, col - 1], axis=0))  #

            for name in unique_parent_taxa:

                if name in color_list:
                    continue

                if name == "NA":
                    color_list[name] = "rgb(0, 0, 0)"

                else:
                    red = np.random.randint(6, high=255)
                    green = np.random.randint(6, high=255)
                    blue = np.random.randint(6, high=255)

                    while red + blue + green > 680:
                        red = np.random.randint(6, high=255)
                        green = np.random.randint(6, high=255)
                        blue = np.random.randint(6, high=255)

                    color_list[name] = "rgb({0}, {1}, {2})".format(
                        red, green, blue)

        if filter_value > 0:
            logging.info("Creating .txt file with filtered results ")
            out_handle = open(filter_title, "w")
            # Counts how many results are filtered
            n_record_filter = 0

        for update_counter, organism in enumerate(taxonomy_counter.keys()):
            ETILI.update_progress(update_counter,
                                  len(taxonomy_counter.keys()))  # Bar progress
            pos_x += 1

            # When there is a filter value > 0 and the counter value is less or the same we will write it in the file
            # and skip the tracing process
            if taxonomy_counter[organism] <= filter_value >= 0:
                logging.info("Writing a new one record in filter file")
                out_handle.write("{0} ==> {1} records\n".format(
                    organism, taxonomy_counter[organism]))
                n_record_filter += 1
                continue

            # If it's the first column we can't group our points based on their parent, because there is no parent
            if col == 0:

                group_taxa = ""
                parent_taxa_text = ""

                if organism == "NA":
                    point_color = "rgb(0, 0, 0)"

                else:
                    red = np.random.randint(6, high=255)
                    green = np.random.randint(6, high=255)
                    blue = np.random.randint(6, high=255)

                    while red + blue + green > 680:  # In order to skip white or extremely light colors
                        red = np.random.randint(6, high=255)
                        green = np.random.randint(6, high=255)
                        blue = np.random.randint(6, high=255)

                    point_color = "rgb({0}, {1}, {2})".format(red, green, blue)

            else:

                if organism == "NA":
                    point_color = "rgb(0, 0, 0)"
                    parent_taxa_text = "({0})".format(
                        ", ".join(unique_parent_taxa_na))
                    group_taxa = "NA"

                else:
                    group_taxa = taxonomy_array[
                        list(taxonomy_array[:, col]).index(organism), col - 1]
                    parent_taxa_text = "({0})".format(group_taxa)
                    point_color = color_list[group_taxa]

            legend_name = '{0} {1}'.format(organism, parent_taxa_text)

            if len(legend_name) > 80:
                legend_name = "{0}...".format(legend_name[0:80])

            if taxonomy_counter[organism] < 4:
                dim = 10

            else:
                dim = 10 + int(taxonomy_counter[organism] * 0.5)

            if dim > 150:
                dim = 150

            if col == 3:
                visi = True
            else:
                visi = False

            logging.debug("Tracing {0} of {1}".format(
                update_counter, len(taxonomy_counter.keys())))
            trace.append(
                go.
                Scattergl(  # Plotly function to trace the point in our scatter plot
                    x=[pos_x],  # A numeric value that increase for every cycle
                    y=[taxonomy_counter[organism]
                       ],  # Counter value of the organism
                    mode=
                    'markers',  # Plot mode, markers means we will trace a point
                    opacity=0.7,
                    text=parent_taxa_text.replace("(", "").replace(")", ""),
                    name=
                    legend_name,  # Value name will be the organism name and the parent in ()
                    visible=visi,
                    hoverinfo="y+text+name",
                    legendgroup=
                    group_taxa,  # This variable will group points based on their parent taxa
                    marker=dict(size=dim,
                                color=point_color,
                                line=dict(color='rgb(255,255,255)', width=1))))
        ETILI.update_progress(len(taxonomy_counter.keys()),
                              len(taxonomy_counter.keys()))

    if filter_value > 0:
        out_handle.close()
        print(
            "\nThere are %i species excluded by your filter, check them in the text file"
            % n_record_filter)
        logging.info(
            "There are %i species with just one entry in the entire database!"
            % n_record_filter)

    print("Preparing to plot...")

    button_steps = []
    # Labels for the steps in the slider
    button_labels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    counter_step = 0

    for i in range(0, len(button_counter)):

        step = dict(
            method='restyle',
            args=['visible', [False] * len(trace)],
            label=button_labels[i],
        )

        for x in range(0, button_counter[i]):
            # Enable all the scatters that are needed each step
            step['args'][1][x + counter_step] = True

        # Add step to step list
        button_steps.append(step)

        counter_step += button_counter[i]

    buttons = [
        dict(
            visible=True,
            type="buttons",
            active=3,
            # font = {'color': 'rgb(107, 255, 109)'},
            # bgcolor='rgb(184, 255, 225)',
            direction='right',
            x=0.5,
            y=-0.1,
            xanchor='center',
            yanchor='bottom',
            buttons=button_steps)
    ]

    layout = go.Layout(title=dict(
        text="{0}<br>(Click legend to toggle traces)".format(title_graph),
        font={'size': 20}),
                       autosize=True,
                       yaxis=dict(type='log', title='N° of sequence'),
                       legend=dict(traceorder='grouped'),
                       hoverlabel=dict(font={'color': 'white'},
                                       namelength=-1,
                                       bordercolor='black'),
                       hovermode='closest',
                       updatemenus=buttons)
    figure = dict(data=trace, layout=layout)

    logging.info(" Plotting the traces")
    plot(figure, filename=plot_file_name)

    print('Choose one of the following options: \n ',
          '1: Plot another file or another query \n ',
          '2: Back to statistical menu \n ', '3: Back to main menu \n ',
          '4: Close the best program you ever used \n ')
    plot_menu_choice = 0

    while plot_menu_choice not in (1, 2, 3, 4):
        try:
            plot_menu_choice = int(input(">> Enter your choice: "))

            if plot_menu_choice == 1:
                ETILI.clear()
                scatterplot()

            elif plot_menu_choice == 2:
                ETILI.clear()
                statistical_module()

            elif plot_menu_choice == 3:
                ETILI.clear()
                ETILI.main_menu()

            elif plot_menu_choice == 4:
                return

            else:
                print("Wrong choice, try again")

        except ValueError:
            print("Error, insert only the numeric value of your choice")


def top10_graph(title, x_values, y_values):
    #y
    yval = y_values
    #x
    xval = x_values

    #plot
    fig = go.Figure([
        go.Bar(x=xval, y=yval, text=xval, textposition='auto', orientation='h')
    ])

    #customize aspect
    fig.update_traces(
        overwrite=True,
        marker={
            "color": xval,
            "colorscale": 'blugrn'
        },
    )

    fig.update_layout(
        title={
            'text': "Top10 " + title + ": ",
            'y': 0.95,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        autosize=False,
        width=900,
        height=950,
        #yaxis_title="",
        font=dict(family="Bodoni 72 Smallcaps", size=17, color='black'),
        paper_bgcolor='#ffffff',
        plot_bgcolor='#e6f0f5',
        bargap=0.4)

    return fig


def worldmap_plot():
    print("--- WORLD MAP MENU ---\n")

    # getting the history of the calls for efetch_call
    efetch_call_settings = ETIDB.database_module("world", 0, None, None,
                                                 [None, None])

    batch_size = 200  # Batch_size value limit our queries to NCBI so we don't get blacklisted
    missing_part = 0
    wm_all = []
    coordinates = [
    ]  # list of coordinates with their respective gene and organism
    countries = [
    ]  # list of countries whenever we don't find coordinates from the record
    no_geo = 0  # Counting how many records didn't put the coordinates nor the country value
    genes = []  # List of genes for the menu on the world map

    if len(efetch_call_settings['query'].replace("(", "").replace(")",
                                                                  "")) < 15:
        print("\nThis is the automatic title: %s" %
              efetch_call_settings['query'].replace("(", "").replace(")", ""))

        try:
            title_choose = int(
                input(
                    "Do you want to keep it? (0 > yes, anything else > no) "))

        except ValueError:
            title_choose = 1
    else:
        title_choose = 1

    if title_choose != 0:
        title_map = input(">> Enter the title of the map: ")

    else:
        title_map = efetch_call_settings['query'].replace("(",
                                                          "").replace(")", "")

    logging.info(' Downloading info for world map')

    print("\nDownloading info from NCBI for world map")
    for start in range(0, efetch_call_settings['counter_id'], batch_size):

        logging.info(' Downloading data for world map: %i of %i' %
                     (start + 1, efetch_call_settings['counter_id']))

        data = ETILRET.efetch_call(start, "nuccore",
                                   efetch_call_settings['counter_id'],
                                   efetch_call_settings['webenv'],
                                   efetch_call_settings['query_key'], "gbc")

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
                INSDQualifiers = INSDSeq.findall(
                    './/INSDSeq_feature-table/INSDFeature/INSDFeature_quals/INSDQualifier'
                )
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
                        while (True):
                            try:
                                if "S" in lat_lon:
                                    lat = 0 - float(
                                        lat_lon[0:lat_lon.index("S") - 1])

                                    if "W" in lat_lon:
                                        lon = 0 - float(
                                            lat_lon[lat_lon.index("S") +
                                                    1:lat_lon.index("W") - 1])
                                        break
                                    else:
                                        lon = float(
                                            lat_lon[lat_lon.index("S") +
                                                    1:lat_lon.index("E") - 1])
                                        break
                                else:
                                    lat = float(lat_lon[0:lat_lon.index("N") -
                                                        1])

                                    if "W" in lat_lon:
                                        lon = 0 - float(
                                            lat_lon[lat_lon.index("N") +
                                                    1:lat_lon.index("W") - 1])
                                        break
                                    else:
                                        lon = float(
                                            lat_lon[lat_lon.index("N") +
                                                    1:lat_lon.index("E") - 1])
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

                wm_all.append({
                    'accession': accession,
                    'org': organism,
                    'country': country,
                    'lat': lat,
                    'lon': lon,
                    'gene': gene
                })

                if lat_lon is not None:
                    coordinates.append({
                        'org': organism,
                        'lat': lat,
                        'lon': lon,
                        'gene': gene
                    })
                    genes.append(gene)

                elif country != "NA":
                    countries.append({
                        'org': organism,
                        'country': country,
                        'gene': gene
                    })

                else:
                    no_geo += 1

    df_wm = pd.DataFrame(wm_all)
    df_wm.to_csv(os.path.join(efetch_call_settings['folder_path'],
                              title_map + '_worldmap.tsv'),
                 sep='\t',
                 index=False)

    #megae dataframe coordinates for plot:
    if coordinates:
        df_coordinates = pd.DataFrame(coordinates)
        #df_coordinates.to_csv(os.path.join(efetch_call_settings['folder_path'],title_map+'_coordinate.tsv'), sep='\t',index=False)
        df_coordinates['count'] = 1
        df_coordinates['gene'] = [
            ','.join(map(str, l)) if type(l) == list else l
            for l in df_coordinates['gene']
        ]
        df_coordinates = df_coordinates.groupby(['lat', 'lon']).agg({
            'org':
            ','.join,
            'gene':
            ','.join,
            'count':
            'sum'
        }).reset_index()

        df_coordinates['gene'] = [
            list(l.split(',')) for l in df_coordinates['gene']
        ]
        df_coordinates['gene'] = [
            dict(Counter(l)) for l in df_coordinates['gene']
        ]
        df_coordinates['gene'] = [
            str(dic).replace('{', '').replace('}', '')
            if type(dic) == dict else dic for dic in df_coordinates['gene']
        ]

        df_coordinates['org'] = [
            list(l.split(',')) for l in df_coordinates['org']
        ]
        df_coordinates['org'] = [
            dict(Counter(l)) for l in df_coordinates['org']
        ]
        df_coordinates['org'] = [
            str(dic).replace('{', '').replace('}', '')
            if type(dic) == dict else dic for dic in df_coordinates['org']
        ]

    #menage dataframe country for plot:
    if countries:
        df_country = pd.DataFrame(countries)
        #df_country.to_csv(os.path.join(efetch_call_settings['folder_path'],title_map+'_paesi.tsv'), sep='\t',index=False)
        df_country['count'] = 1
        df_country['gene'] = [
            ','.join(map(str, l)) if type(l) == list else l
            for l in df_country['gene']
        ]
        df_country = df_country.groupby(['country']).agg({
            'org': ','.join,
            'gene': ','.join,
            'count': 'sum'
        }).reset_index()

        df_country['gene'] = [list(l.split(',')) for l in df_country['gene']]
        df_country['gene'] = [dict(Counter(l)) for l in df_country['gene']]
        df_country['gene'] = [
            str(dic).replace('{', '').replace('}', '')
            if type(dic) == dict else dic for dic in df_country['gene']
        ]

        df_country['org'] = [list(l.split(',')) for l in df_country['org']]
        df_country['org'] = [dict(Counter(l)) for l in df_country['org']]
        df_country['org'] = [
            str(dic).replace('{', '').replace('}', '')
            if type(dic) == dict else dic for dic in df_country['org']
        ]
        #df_country.to_csv(os.path.join(efetch_call_settings['folder_path'],title_map+'onlycountry_worldmap.tsv'), sep='\t',index=False)

    ETILI.update_progress(efetch_call_settings['counter_id'],
                          efetch_call_settings['counter_id'])
    cities = []
    names = []
    element_count = 0
    print("")
    print(
        "Of %i records on NCBI, %i has no info about the location of sampling."
        % (efetch_call_settings['counter_id'], no_geo))
    time.sleep(2)
    print("\nWorld map loading...")

    if len(coordinates) > 0:
        for element in df_coordinates.itertuples():
            element_count += 1
            #update_progress(element_count, len(df_coordinates))

            if element.org not in names:
                names.append(element.org)
                showlegend = True
            else:
                showlegend = False

            textspecies = 'Accession count :'
            if len(list(element.org.split(','))) > 5:
                chunks = [
                    list(element.org.split(','))[x:x + 5]
                    for x in range(0, len(list(element.org.split(','))), 5)
                ]
                print(chunks)
                for chunk in chunks:
                    textspecies += str(chunk)
                    textspecies += '<br>'
            else:
                textspecies += str([element.org])

            textgene = 'Gene count :'
            if len(list(element.gene.split(','))) > 5:
                chunks = [
                    list(element.gene.split(','))[x:x + 5]
                    for x in range(0, len(list(element.gene.split(','))), 5)
                ]
                print(chunks)
                for chunk in chunks:
                    textgene += str(chunk)
                    textgene += '<br>'
            else:
                textgene += str([element.gene])

            namespecies = ''
            if len(list(element.org.split(','))) > 3:
                chunks = [
                    list(element.org.split(','))[x:x + 3]
                    for x in range(0, len(list(element.org.split(','))), 3)
                ]
                print(chunks)
                for chunk in chunks:
                    namespecies += str(chunk)
                    namespecies += '<br>'
            else:
                namespecies += str([element.org])

            cities.append(
                go.Scattergeo(
                    locationmode='country names',
                    lon=[element.lon],
                    lat=[element.lat],
                    text=textspecies + '<br>' + textgene,
                    #legendgroup=element['gene'],
                    showlegend=showlegend,
                    marker=go.scattergeo.Marker(
                        symbol=4,
                        opacity=0.7,
                        size=10,
                        color='rgb(88,214,141)',
                        line=go.scattergeo.marker.Line(  # Contorno punto
                            width=1, color='rgb(40,40,40)'),
                        sizemode='area'),
                    name=namespecies))

    if len(countries) > 0:

        for element in df_country.itertuples():
            element_count += 1
            #update_progress(element_count, len(df_country))

            if element.count not in names:
                names.append(element)
                showlegend = True
            else:
                showlegend = False

            textspecies = 'Accession count :'
            if len(list(element.org.split(','))) > 5:
                chunks = [
                    list(element.org.split(','))[x:x + 5]
                    for x in range(0, len(list(element.org.split(','))), 5)
                ]
                print(chunks)
                for chunk in chunks:
                    textspecies += str(chunk)
                    textspecies += '<br>'
            else:
                textspecies += str([element.org])

            textgene = 'Gene count :'
            if len(list(element.gene.split(','))) > 5:
                chunks = [
                    list(element.gene.split(','))[x:x + 5]
                    for x in range(0, len(list(element.gene.split(','))), 5)
                ]
                print(chunks)
                for chunk in chunks:
                    textgene += str(chunk)
                    textgene += '<br>'
            else:
                textgene += str([element.gene])

            cities.append(
                go.Scattergeo(
                    locationmode='country names',
                    locations=[element.country],
                    text=textspecies + '<br>' + textgene,
                    # legendgroup=element['gene'],
                    showlegend=showlegend,
                    marker=go.scattergeo.Marker(
                        opacity=0.7,
                        symbol=0,
                        size=10,
                        color='rgb(203,67,53)',
                        line=go.scattergeo.marker.Line(  # Contorno punto
                            width=1, color='rgb(40,40,40)'),
                        sizemode='area'),
                    name=str('Total accession: ') + str(element.count) +
                    str('-') + str(element.country)))


###modify layout extaxsi:
    layout = go.Layout(
        title=go.layout.Title(text="{0}{1}".format(
            title_map, '<br>(Click legend to toggle traces)')),
        showlegend=True,
        hoverlabel=dict(  # font={'color': 'white'},
            namelength=-1,
            # bordercolor='black'
        ),
        # legend=dict(traceorder='grouped'),
        geo=go.layout.Geo(
            scope='world',
            projection=go.layout.geo.Projection(type='natural earth'),
            showland=True,
            landcolor='#DEC86B',
            oceancolor='#69A8FF',
            showocean=True,
            showrivers=True,
            showlakes=True,
            showcountries=True,
            subunitwidth=1,
            countrywidth=1,
            subunitcolor="rgb(255, 255, 255)",
            countrycolor="rgb(255, 255, 255)"))

    fig = go.Figure(data=cities, layout=layout)
    plot(fig,
         filename=ETIL.rename_file(efetch_call_settings['folder_path'],
                                   title_map, "_worldmap.html"))

    if missing_part == 0:
        print("\n ---- World map plot has been created. ----\n")
        time.sleep(30)

    else:
        print(
            "\n ---- World map plot has been created but some info are due to NCBI server errors. ----"
        )
        print(" ---- Check log file for more details. ---- \n")


def sunburst_plot():
    print("Choose which type of input you will use\n"
          "Enter one of the following options number\n"
          "  1 - Taxonomy file from this tool\n"
          "  2 - Doing a query on NCBI and directly plot it\n"
          "  3 - Manual created file")

    while True:
        try:
            query_or_path = int(
                input("\n>>Enter your choice (only the number): "))

            while query_or_path not in (1, 2, 3):
                print(ETILI.color.RED + "Enter a valid choice, try again. \n" +
                      ETILI.color.END)
                query_or_path = input(
                    ">> Enter your choice (only the number): ")

            break

        except ValueError:
            print(ETILI.color.RED + "Please enter only the numeric value" +
                  ETILI.color.END)
            print("Please enter only the numeric value")

    # Until they put an existing file path, it will continue to ask
    while True:

        if query_or_path == 2:
            plot_file_path = ETIDB.database_module("scatter", 0, None, None, [
                None, None
            ])  # A call to download the taxonomy file and returns the path

        else:
            logging.info(" Getting the taxonomy file path ")
            plot_file_path = str(
                input(">> Enter the taxonomy file path name : "))

        try:
            if plot_file_path[-3:] == 'csv':
                with open(plot_file_path, 'r') as axe:
                    r = csv.reader(axe, delimiter=",")
                    r_list = list(r)
                    break

            elif plot_file_path[-3:] == 'tsv':
                with open(plot_file_path, 'r') as axe:
                    r = csv.reader(axe, delimiter="\t")
                    r_list = list(r)
                    break
            else:
                print(
                    "File extension not CSV or TSV, please use the right format."
                )

        except IOError:
            logging.warning(
                " No file found or permission denied, path entered: %s" %
                plot_file_path)
            print(
                ETILI.color.RED +
                "No file found or permission denied, please check your file or location"
                + ETILI.color.END)
            print("File location: ", plot_file_path, "\n")

    taxonomy_list = []
    taxa_level = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    print(r_list[0])
    if input("Is this row a title? (Y > yes, N > no) ") in ("Y", "y"):
        skip_first_row = True
        print("Ok, i'll skip it")

    else:
        print("Good!\n")
        skip_first_row = False

    for row in r_list:

        if skip_first_row:  # Excluding the first row when we find titles
            skip_first_row = False
            continue

        if query_or_path == 3:
            taxonomy_list.append(row)
        else:
            taxonomy_list.append(row[1].split(";"))

    taxonomy_array = np.array(taxonomy_list)

    print(
        "Now enter the minimum results value needed to be plotted for each taxonomy"
    )
    filter_value = None
    while filter_value is None:

        try:
            filter_value = int(input('>> Enter filter value: '))
            if filter_value < 0:
                print(ETILI.color.RED + "Enter a valid positive number!" +
                      ETILI.color.END)
                filter_value = None

        except ValueError:
            print(ETILI.color.RED + "Enter a number!" + ETILI.color.END)

    title_graph = input(">> Enter the title for the plot: ")  #

    while "/" in plot_file_path:
        x = plot_file_path.find("/")
        plot_file_path = plot_file_path[x + 1:]

    if "." in plot_file_path:
        x = plot_file_path.find(".")
        plot_file_path = plot_file_path[:x]

    plot_file_name = ETIL.rename_file(directory, title_graph, '.html')
    filter_title = ETIL.rename_file(
        directory,
        "{0}_with_less_than_{1}_results".format(plot_file_path,
                                                str(filter_value)), ".txt")
    if os.path.exists(plot_file_name):
        x = 0
        overwrite = int(
            input(
                ETILI.color.YELLOW +
                "Exist already a file with the same name, do you want to overwrite it? (0 = yes, 1 = no)"
                + ETILI.color.END))
    else:
        overwrite = 0

    while os.path.exists(plot_file_name) and overwrite != 0:
        x += 1
        plot_file_name = ETIL.rename_file(directory, title_graph,
                                          '_%i.html' % x)
        filter_title = ETIL.rename_file(
            directory,
            "{0}_with_less_than_{1}_results".format(title_graph,
                                                    str(filter_value)),
            "_%i.txt" % x)

    while True:
        try:
            depth = int(
                input(
                    "Smaller depths = light plots, default level is 2, available values: 2, 3, 4, 5, 6"
                    + "\nEnter the depth of your sunburst : "))
            if depth not in (2, 3, 4, 5, 6):
                print(
                    "Please enter an natural number not less than 2 not more than 6!"
                )
            break

        except ValueError:
            print("Please insert a numeric value from 2 to 6.")

    ##############
    # Initializing the list of traces for sunburst
    lab = [" "]  # List of labels
    par = [""]  # List of parents
    val = [0]  # List of their values

    lab.append("NA")
    par.append(" ")
    val.append(0)

    for col in range(len(taxonomy_list[0])):

        print("\nElaborating {0}{1}{2}".format(ETILI.color.BOLD,
                                               taxa_level[col],
                                               ETILI.color.END))
        logging.info("Elaborating data of column: %i" % int(col))
        taxonomy_counter = collections.Counter(taxonomy_array[:, col])

        if filter_value > 0:
            logging.info("Creating .txt file with filtered results ")
            out_handle = open(filter_title, "w")
            n_record_filter = 0  # Counts how many results are filtered

        for update_counter, organism in enumerate(taxonomy_counter.keys(
        )):  # With enumerate i can track the work done
            ETILI.update_progress(update_counter,
                                  len(taxonomy_counter.keys()))  # Bar progress

            # When there is a filter value > 0 and the counter value is less or the same we will write it in the file
            # and skip the tracing process
            if taxonomy_counter[organism] <= filter_value >= 0:
                logging.info("Writing a new one record in filter file")
                out_handle.write("{0} ==> {1} records\n".format(
                    organism, taxonomy_counter[organism]))
                n_record_filter += 1
                continue

            if organism == "NA":
                val[1] += taxonomy_counter[organism]
                continue

            lab.append(organism)
            val.append(taxonomy_counter[organism])

            # If it's the first column we can't group our points based on their parent, because there is no parent
            if col == 0:
                par.append(" ")

            else:
                x = 1
                while col >= x:
                    if taxonomy_array[
                            list(taxonomy_array[:, col]).index(organism),
                            col - x] == "NA":
                        x += 1

                    else:
                        break

                if col >= x:
                    par.append(taxonomy_array[
                        list(taxonomy_array[:, col]).index(organism), col - x])

                else:
                    par.append(" ")

        ETILI.update_progress(len(taxonomy_counter.keys()),
                              len(taxonomy_counter.keys()))

    if val[1] == 0:
        lab.pop(1)
        par.pop(1)
        val.pop(1)

    fig = go.Figure(
        go.Sunburst(
            labels=lab,
            parents=par,
            values=val,
            # branchvalues="total",
            hovertemplate='<b>%{label}</b><br>N° records: %{value}<br>'
            +  # %{percentEntry:.2%} of the parent<br> strange values
            '<extra></extra>',
            maxdepth=depth,
        ))

    if filter_value > 0:
        out_handle.close()
        print(
            "\nThere are %i species excluded by your filter, check them in the text file"
            % n_record_filter)
        logging.info(
            "There are %i species with just one entry in the entire database!"
            % n_record_filter)
        time.sleep(2)

    logging.info(" Plotting the traces")
    fig.update_layout(margin=dict(t=0, l=0, r=0, b=0))
    plot(fig, filename=plot_file_name)

    print('\nChoose one of the following options: \n ',
          '1: Another plot from file or query \n ',
          '2: Back to statistical menu \n ', '3: Back to main menu \n ',
          'CTRL + C: Close the best program you ever used \n ')
    plot_menu_choice = 0

    while plot_menu_choice not in (1, 2, 3, 4):
        try:
            plot_menu_choice = int(input(">> Enter your choice: "))

            if plot_menu_choice == 1:
                ETILI.clear()
                sunburst_plot()

            elif plot_menu_choice == 2:
                ETILI.clear()
                statistical_module()

            elif plot_menu_choice == 3:
                ETILI.clear()
                ETILI.main_menu()

            elif plot_menu_choice == 4:
                return

            else:
                print("Wrong choice, try again")

        except ValueError:
            print("Error, insert only the numeric value of your choice")  # #
