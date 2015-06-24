#!/usr/bin/env python

from distutils.core import setup

setup(name='ortho',
      version='0.5',
      description='miscellany of python packages for all and sundry analyses',
      author='Peter Oxley',
      author_email='oxpeter+git@gmail.com',
      url='https://oxpeter@bitbucket.org/oxpeter/ortho.git',
      package_data={'genomepy': ['data/*.cfg'], 'genomepy': ['data/README']},
      py_modules=[  'blastfaster', 'brain_machine', 'degrees',
                    'kegg', 'longestORF',
                    'plate_analysis', 'textmod',
                  ],
      requires=[ 'argparse' ]
     )