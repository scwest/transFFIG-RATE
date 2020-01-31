from setuptools import setup

import subprocess
import os

setup(
      name = 'transFIGG-RATE',
      version = '0.1',
      description = 'Generation of intermediate (between gene and transcript) granules.',
      author = 'Sean West',
      url = 'tbd',
      packages = ['transfigg_prep', 'transfigg_rate'],
      package_data = {
                      '':['docs/*']
                      },
      entry_points = {
                      'console_scripts':[
                                         'transfigg_rate = transfigg_rate.control:smain',
                                         'transfigg_prep = transfigg_prep.control:smain',
                                         'transfigg_clustalo = transfigg_prep/transfigg_clustalo.py',
                                         'transfigg_mafft = transfigg_prep/transfigg_mafft.py',
                                         'transfigg_muscle = transfigg_prep/transfigg_muscle.py',
                                         'transfigg_tcoffee = transfigg_prep/transfigg_tcoffee.py',
                                         
                                         ]
                      },
      long_description = 'Generation of intermediate (between gene and transcript) granules.',
      classifiers = ['Programming Language :: Python', \
                     'Programming Language :: Python :: 3', \
                     'Operating System :: Unix', \
                     'Development Status :: 1 - Planning', \
                     'Intended Audience :: Science/Research', \
                     'Topic :: Scientific/Engineering :: Bio-Informatics', \
                     'Topic :: Scientific/Engineering :: Mathematics', \
                     'License :: OSI Approved :: GNU General Public License v3 (GPLv3)', \
                     'Natural Language :: English']
      )










