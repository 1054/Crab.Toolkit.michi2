#!/usr/bin/env python
# 
# This code will subtract spectral line flux contribution to the SED photometry flux densities. 
# Inputs: redshift, 
# 
# 

import os, sys, re
import numpy as np
import astropy
from astropy.table import Table
from copy import copy



#####################################
#               USAGE               #
#####################################

def usage():
    print('Usage: ')
    print('    michi2_subtract_spectral_line_flux_contribution.py datatable_photometry.txt \\')
    print('                                                       -redshift 3.0 \\')
    print('                                                       -filter-name "^ALMA.*" \\')
    print('                                                       -filter-name-column 5 \\')
    print('                                                       -freq-support 667.84 669.84 669.61 671.61 671.42 673.42 673.22 675.22 \\')
    print('                                                       -line-name CO CI NII CII H2O \\')
    print('                                                       -IR-luminosity 5e12')
    print('                                                       [-known-line-name-and-flux \"CO(2-1)\" \"3.0 Jy km/s\"]')
    print('                                                       [-known-line-name-and-flux \"CO(3-2)\" \"5.0 Jy km/s\"]')
    print('Notes:')
    print('    The input file name "datatable_photometry.txt" must be put in front.')
    print('    The input table "datatable_photometry.txt" must have wavelength (um) in the first column.')
    print('    The filter name is an option to subtract line contamination for only selected filters. The input requires Python regex string.')
    print('    The filter name column indicates which column is the filter name in the input table "datatable_photometry.txt", which is the 5th column in default.')



####################################
#               MAIN               #
####################################

if __name__ == '__main__':
    
    # Read user input
    input_table_file = ''
    input_redshift = np.nan
    input_freq_support = []
    input_IR_luminosity = np.nan
    input_line_names = []
    input_filter_names = []
    input_filter_regexes = []
    input_filter_name_column = 5
    speed_of_light_kms = 2.99792458e5
    iarg = 0
    tstr = ''
    tmode = ''
    isopt = False
    while iarg < len(sys.argv):
        tstr = re.sub(r'^[-]+', r'', sys.argv.lower()) # lower case current input argument string
        isopt = False # whether current input argument is an option (starting with "-")
        if tstr == 'redshift' or tstr == 'z':
            tmode = 'redshift'
            isopt = True
        elif tstr == 'IR-luminosity':
            tmode = 'IR_luminosity'
            isopt = True
        elif tstr == 'freq-support' or tstr == 'freq':
            tmode = 'freq_support'
            isopt = True
        elif tstr == 'line-names' or tstr == 'line-name' or tstr == 'linenames' or tstr == 'linename':
            tmode = 'line_names'
            isopt = True
        elif tstr == 'filter-names' or tstr == 'filter-name' or tstr == 'filternames' or tstr == 'filtername':
            tmode = 'filter_names'
            isopt = True
        elif tstr == 'filter-name-column':
            tmode = 'filter_name_column'
            isopt = True
        else:
            if tmode == '':
                input_table_file = sys.argv[iarg]
        # 
        if not isopt:
            if tmode == 'redshift':
                input_redshift = float(sys.argv[iarg])
            elif tmode == 'IR_luminosity':
                input_IR_luminosity = float(sys.argv[iarg])
            elif tmode == 'freq_support':
                input_freq_support.append(float(sys.argv[iarg]))
            elif tmode == 'line_names':
                input_line_names.append(sys.argv[iarg])
            elif tmode == 'filter_names':
                input_filter_names.append(sys.argv[iarg])
                input_filter_regexes.append(re.compile(sys.argv[iarg]))
        # 
        iarg = iarg + 1
    
    # Check user input
    if input_table_file == '' or np.isnan(input_redshift) or len(input_freq_support) == 0 or np.isnan(input_IR_luminosity):
        usage()
        sys.exit()
    
    # Read data table
    tb = Table.read(input_table_file, format='ascii')
    
    # store freq support as np.array
    freq_support = np.array(input_freq_support)
    
    # Find out the intersected frequencies/wavelengths
    wavelengths = tb.field(tb.colnames[0]).data # um
    frequencies = speed_of_light_kms / wavelengths # GHz
    filternames = [x for x in tb['FilterName'].data]
    for i in range(len(frequencies)):
        
        if tb['FilterName'].startswith('ALMA'):
            if freq_support.intersect1d():
                line_freqs, line_names = find_radio_lines_in_frequency_range(input_freq_support, Redshift=input_redshift, include_faint_lines = False)
        







