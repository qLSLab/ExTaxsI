from distutils.core import setup
setup(
  name = 'extaxsi',
  packages = ['extaxsi'],
  version = '0.1',
  license='MIT',
  description = 'Extaxsi is a bioinformatic library aimed to elaborate and visualize molecular and taxonomic informations',
  author = 'alberto brusati',
  author_email = 'alberto.brusati@gmail.com',
  url = 'https://github.com/user/extaxsi',
  download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['bioinformatic', 'ncbi', 'molecular data', 'visualization', 'taxonomy', 'converter'],
  install_requires=[            # I get to this in a second
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
          'psutil',
          'requests'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
  ],
)
