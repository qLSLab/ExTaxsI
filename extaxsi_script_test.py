import extaxsi

extaxsi.load_configurations("user@email.com","api_key_here", taxa_database_update = 'no')

extaxsi.db_creation(text_search='txid8832',
            accession_taxonomy_output=True,
            fasta_output=False,
            marker_output=True,
            top10_plot=True,
            enrich_output=True)

extaxsi.db_creation(file_search='examples/A_accession_list_example.tsv',
            input_file_type='A',
            accession_taxonomy_output=True,
            fasta_output=True,
            marker_output=True,
            top10_plot=True)

extaxsi.taxonomyID_converter(file_search = 'examples/O_organism_list_example.tsv', input_type = 'O')

extaxsi.sunburst_plot("download/txid8832_taxonomy.tsv", "example")

extaxsi.scatterplot("download/txid8832_taxonomy.tsv", "example")

extaxsi.worldmap_plot("download/txid8832_enriched.tsv",'example_worldmap')
