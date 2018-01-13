#!/usr/bin/env python3.6
# 

import os
import sys
import numpy
import astropy
import astropy.io.ascii as asciitable
from copy import copy




####################################
#               MAIN               #
####################################

if not len(sys.argv) > 2:
    print('Usage: michi2_filter_flux_3sigma.py input_flux.txt output_flux.txt')
    sys.exit()

data_table = asciitable.read(sys.argv[1])

if not len(data_table.colnames) >= 3:
    print('Error! The input flux data table does not have at least three columns: wavelength, flux density and error in flux density.')
    sys.exit()

w = data_table.field(data_table.colnames[0])
f = data_table.field(data_table.colnames[1])
ferr = data_table.field(data_table.colnames[2])
mask = (f<2.0*ferr)
isel = numpy.argwhere(mask).flatten()
print(isel)
print(data_table)
data_table.remove_rows(isel)
print(data_table)

out_file = sys.argv[2]
asciitable.write(data_table, out_file, Writer=asciitable.FixedWidthTwoLine, overwrite=True)
os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(out_file))
os.system('sed -i.bak -e "2d" "%s"'%(out_file))









