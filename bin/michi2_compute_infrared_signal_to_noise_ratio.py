#!/usr/bin/env python
# 

import os
import sys
import numpy
import astropy
import astropy.io.ascii as asciitable
from copy import copy
import json




####################################
#               MAIN               #
####################################

if not len(sys.argv) > 2:
    print('Usage: michi2_compute_infrared_signal_to_noise_ratio.py input_flux.txt output_file.txt')
    sys.exit()

data_table = asciitable.read(sys.argv[1])

if not len(data_table.colnames) >= 3:
    print('Error! The input flux data table does not have at least three columns: wavelength, flux density and error in flux density.')
    sys.exit()

# select infrared data points
w = data_table.field(data_table.colnames[0])
f = data_table.field(data_table.colnames[1])
ferr = data_table.field(data_table.colnames[2])

# set zero error data error to 1/3 flux
mask = (ferr==0)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    ferr[mask] = f[mask] / 3.0
    for iseli in isel:
        print('Limited row %d error from zero to 1/3: w = %s, f = %s, ferr = %s'%(iseli, w[iseli], f[iseli], ferr[iseli]))

# compute S/N
snr = {}
snr['SNR_8_1000_um'] = 0.0
snr['SNR_100_1000_um'] = 0.0
snr['SNR_100_3000_um'] = 0.0

mask = (w>=8.0) & (w<=1000.0)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    snr_array = (f[mask]/ferr[mask])**2
    snr['SNR_8_1000_um'] = numpy.sqrt(numpy.sum(snr_array))
    snr['NUM_8_1000_um'] = len(snr_array)

mask = (w>=100.0) & (w<=1000.0)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    snr_array = (f[mask]/ferr[mask])**2
    snr['SNR_100_1000_um'] = numpy.sqrt(numpy.sum(snr_array))
    snr['NUM_100_1000_um'] = len(snr_array)

mask = (w>=100.0) & (w<=3000.0)
isel = numpy.argwhere(mask).flatten()
if len(isel) > 0:
    snr_array = (f[mask]/ferr[mask])**2
    snr['SNR_100_3000_um'] = numpy.sqrt(numpy.sum(snr_array))
    snr['NUM_100_3000_um'] = len(snr_array)

# output
out_file = sys.argv[2]
with open(out_file, 'w') as fp:
    json.dump(snr, fp)
    fp.write('\n')
    fp.close()

os.system('cat "%s"'%(out_file))
print('Output to "%s"!'%(out_file))








