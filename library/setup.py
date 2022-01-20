from distutils.core import setup

setup(
  name = 'extaxsi',
  packages = ['extaxsi'],
  version = '0.3.9',
  license='MIT',
  description = 'Extaxsi is a bioinformatic library aimed to elaborate and visualize molecular and taxonomic informations',
  long_description="""ExTaxsI functions \n

Library overview \n

ExTaxsI functions can be used separately from the ExTaxsI tool. In general, there are 6 main functions for the following tasks: \n


*   load_configurations() - **mandatory before using ExTaxsi functions** - setting NCBI API key to allow data query \n
*   db_creation() - create FASTA, accession with taxonomy or enriched database \n
*   taxonomyID_converter() - convert taxid into organism taxonomy and viceversa \n
*   sunburst_plot() - create sunburst plot \n
*   scatterplot() - create scatterplot \n
*   worldmap_plot() - create worldmap plot \n

In the [github](https://github.com/qLSLab/ExTaxsI/blob/master/library/functions_installation_and_tutorial.ipynb) notebook installation and tutorials for configuration and usage are provided. \n
Input files used are available in the **examples** directory. \n

			For any doubt, please write to g.agostinetto@campus.unimib.it""",
  long_description_content_type="text/markdown",
  author = 'alberto brusati',
  author_email = 'alberto.brusati@gmail.com',
  url = 'https://github.com/user/extaxsi',
  download_url = 'https://github.com/qLSLab/ExTaxsI/archive/refs/tags/v_03.tar.gz',
  keywords = ['bioinformatic', 'ncbi', 'molecular data', 'visualization', 'taxonomy', 'converter'],
  install_requires=[
          'biopython',
          'numpy',
          'scipy',
          'matplotlib',
          'ipython',
          'pandas',
          'sympy',
          'nose',
          'plotly',
          'ete3',
          'kaleido'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
  ],
)
