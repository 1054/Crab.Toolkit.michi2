#!/usr/bin/env python
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

# 
w = data_table.field(data_table.colnames[0])
f = data_table.field(data_table.colnames[1])
ferr = data_table.field(data_table.colnames[2])
mask = (f<2.0*ferr) | (w>10000)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    #print(isel)
    #print(data_table)
    data_table.remove_rows(isel)
    #print(data_table)

# deal with duplicated w
i = 0
while i < len(data_table):
    w = data_table.field(data_table.colnames[0])
    f = data_table.field(data_table.colnames[1])
    ferr = data_table.field(data_table.colnames[2])
    # identify duplicated w
    mask2 = (w==w[i])
    isel2 = numpy.argwhere(mask2).flatten()
    if len(isel2) >= 2:
        # found duplicated w
        print('Found wavelength-duplicated rows: %s'%(isel2))
        print(data_table[isel2])
        f_to_average = f[mask2]
        ferr_to_average = ferr[mask2]
        f_averaged = numpy.sum(f_to_average/ferr_to_average**2)/numpy.sum(1/ferr_to_average**2)
        ferr_averaged = numpy.sqrt(1/numpy.sum(ferr_to_average**(-2))) # error propagation of weighted mean, see -- http://www.physics.umd.edu/courses/Phys261/F06/ErrorPropagation.pdf
        # limit S/N not larger than 10
        #if ferr_averaged < f_averaged/10.0:
        #    ferr_averaged = f_averaged/10.0
        # store into data_table
        f[i] = f_averaged # change f will directly change data_table!
        ferr[i] = ferr_averaged # change ferr will directly change data_table!
        print('Averaged wavelength-duplicated rows: w = %s, f = %s, ferr = %s'%(w[i], f_averaged, ferr_averaged))
        # remove those duplicated rows, but keep current i row.
        isel3 = isel2[(isel2 != i)]
        for iseli in isel3:
            print('data_table.remove_rows(%d)'%(iseli))
        data_table.remove_rows(isel3)
    i = i+1

# limit S/N to be not larger than 10
w = data_table.field(data_table.colnames[0])
f = data_table.field(data_table.colnames[1])
ferr = data_table.field(data_table.colnames[2])
mask = (ferr<f/10.0)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    ferr[mask] = f[mask] / 10.0
    for iseli in isel:
        print('Limited row %d S/N no larger than 10: w = %s, f = %s, ferr = %s'%(iseli, w[iseli], f[iseli], ferr[iseli]))

# output
out_file = sys.argv[2]
asciitable.write(data_table, out_file, Writer=asciitable.FixedWidthTwoLine, overwrite=True)
os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(out_file))
os.system('sed -i.bak -e "2d" "%s"'%(out_file))
if os.path.isfile(out_file+'.bak'):
    os.system('rm "%s"'%(out_file+'.bak'))
print('')








