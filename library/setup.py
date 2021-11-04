from distutils.core import setup
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
  name = 'extaxsi',
  packages = ['extaxsi'],
  version = '0.3.3',
  license='MIT',
  description = 'Extaxsi is a bioinformatic library aimed to elaborate and visualize molecular and taxonomic informations',
  long_description=long_description,
  long_description_content_type="text/markdown",
  author = 'alberto brusati',
  author_email = 'alberto.brusati@gmail.com',
  url = 'https://github.com/user/extaxsi',
  download_url = 'https://github.com/user/extaxsi/archive/refs/tags/v_03.tar.gz',
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
