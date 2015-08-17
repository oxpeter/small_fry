#!/usr/bin/env python

from distutils.core import setup

setup(name='small_fry',
      version='0.5',
      description='miscellany of python packages for all and sundry analyses',
      author='Peter Oxley',
      author_email='oxpeter+git@gmail.com',
      url='https://oxpeter@bitbucket.org/oxpeter/ortho.git',
      package_data={'': ['data/*.txt'], '': ['data/*.gff']},
      py_modules=[  'blastfaster',
                    'brain_machine',
                    'degrees',
                    'kegg',
                    'longestORF',
                    'compare_clusters',
                    'orthotree',
                    'plate_analysis',
                    'textmod',
                  ],
      requires=[ 'argparse' ]
     )