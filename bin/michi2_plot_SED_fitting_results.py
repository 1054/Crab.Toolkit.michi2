#!/usr/bin/env python3.6
# 

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')
sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabplot')

from CrabTable import *
from CrabPlot import *

import glob
import math
import numpy
import astropy
from astropy import units
from astropy.io import fits
import re
import json





#########################################
#               Functions               #
#########################################











##########################################
#               MAIN PROGRAM             #
##########################################

if len(sys.argv) == 1:
    
    print('Usage: michi2_plot_SED_fitting_results.py fit_5.out')
    sys.exit()

else:
    # 
    # Read chi2 table
    DataFile = sys.argv[1]
    print('# Reading "%s"'%(DataFile))
    DataTable = CrabTable(DataFile, verbose=0)
    # 
    # Read fit info
    InfoFile = DataFile + '.info'
    if os.path.isfile(InfoFile):
        print('# Reading "%s"'%(InfoFile))
        InfoDict = CrabTableReadInfo(InfoFile, verbose=0)
    CheckInfoDictOK = True
    for InfoKey in ['OBS', 'NLIB', 'OUT']:
        if InfoKey not in InfoDict:
            print('Error! Key "%s" is not in the InfoFile "%s"!'%(InfoKey, InfoFile))
            CheckInfoDictOK = False
    if CheckInfoDictOK is False:
        sys.exit()
    # 
    print(InfoDict)
    # 
    # Sort chi2 table
    print(DataTable.headers)
    












