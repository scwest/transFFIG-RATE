from setuptools import setup

import subprocess
import os

setup(
      name = 'transFFIG-RATE',
      version = '0.1',
      description = 'Generation of intermediate (between gene and transcript) granules.',
      author = 'Sean West',
      url = 'tbd',
      packages = ['transffig_prep', 'transffig_rate'],
      package_data = {
                      '':['docs/*']
                      },
      entry_points = {
                      'console_scripts':[
                                         'transffig_rate = transffig_rate.control:smain',
                                         'transffig_prep = transffig_prep.control:smain',
                                         'transffig_clustalo = transffig_clustalo.transffig_clustalo:smain',
                                         'transffig_mafft = transffig_mafft.transffig_mafft:smain',
                                         'transffig_muscle = transffig_muscle.transffig_muscle:smain',
                                         'transffig_tcoffee = transffig_tcoffee.transffig_tcoffee:smain'
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










