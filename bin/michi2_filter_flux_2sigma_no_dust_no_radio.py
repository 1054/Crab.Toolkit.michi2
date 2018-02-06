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
mask = (f<2.0*ferr) | (w>20)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    #print(isel)
    #print(data_table)
    data_table.remove_rows(isel)
    #print(data_table)

# set zero error data error to 1/3 flux
w = data_table.field(data_table.colnames[0])
f = data_table.field(data_table.colnames[1])
ferr = data_table.field(data_table.colnames[2])
mask = (ferr==0)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    ferr[mask] = f[mask] / 3.0
    for iseli in isel:
        print('Limited row %d error from zero to 1/3: w = %s, f = %s, ferr = %s'%(iseli, w[iseli], f[iseli], ferr[iseli]))

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

# sort
data_table.sort(data_table.colnames[0])

# output
out_file = sys.argv[2]
asciitable.write(data_table, out_file, Writer=asciitable.Ipac, delimiter='    ', overwrite=True)
#asciitable.write(data_table, sys.stdout, Writer=asciitable.Ipac, delimiter='  ')
with open(out_file, 'r+') as fp:
    out_content = fp.readlines() # read everything in the file
    out_iline = 0
    out_header = [] # Ipac format has multiple comment lines (commented by the char '\\') and 4 header lines.
    fp.seek(0)
    while out_iline < len(out_content):
        if out_content[out_iline][0] == '\\':
            # if his is a commented line, then we change the comment mark to '#'
            out_content[out_iline] = '#' + out_content[out_iline][1:]
            fp.write(out_content[out_iline])
        else:
            if len(out_header) == 0:
                # if this is the first header line, then replace the first white space by '#', or if there is no white space, preprend '#'.
                if out_content[out_iline][0] == ' ':
                    out_content[out_iline] = '#' + out_content[out_iline][1:]
                else:
                    out_content[out_iline] = '#' + out_content[out_iline]
                # append header to 'out_header' list
                out_header.append(out_content[out_iline])
                # write only one header line
                fp.write(out_content[out_iline])
                # 
            elif len(out_header) < 4:
                # append header to 'out_header' list
                out_header.append(out_content[out_iline])
                # skip the 2nd to 4th header line
                pass
            else:
                # write data line
                fp.write(out_content[out_iline])
                # 
        out_iline = out_iline + 1
    fp.truncate()
    fp.close()
#os.system('sed -i.bak -e "$(grep \"\\\" %s | wc -l)s/^ /#/" "%s"'%(out_file, out_file))
#os.system('sed -i.bak -e "2d;3d;4d" "%s"'%(out_file))
#if os.path.isfile(out_file+'.bak'):
#    os.system('rm "%s"'%(out_file+'.bak'))
print('Output to "%s"!'%(out_file))








