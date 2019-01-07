#!/usr/bin/env python
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
    print('    michi2_convert_column_table_to_row_table.py datatable_photometry_100_columns.txt datatable_photometry_100_rows.txt')
    print('Notes:')
    print('    # The first input should be a multi-column data table, with first column being ID, second column being redshift, then next columns ')
    print('      should be flux and error of each photometry.')
    print('    # The second input should be the output multi-row data table with different photometry bands in different rows.')
    print('    # If the input table has more than one row, each row will be output to a differnt output table. ')



####################################
#               MAIN               #
####################################

if __name__ == '__main__':
    
    # Prepare parameters
    input_table_file = ''
    output_table_file = ''
    input_redshift = np.nan
    
    # Read user input
    iarg = 1
    tstr = ''
    while iarg < len(sys.argv):
        tstr = re.sub(r'^[-]+', r'', sys.argv[iarg].lower())
        if tstr == 'redshift' or tstr == 'z':
            input_redshift = float(sys.argv[iarg])
        elif input_table_file == '':
            input_table_file = sys.argv[iarg]
        elif output_table_file == '':
            output_table_file = sys.argv[iarg]
        # 
        iarg = iarg + 1
    
    # Check user input
    if input_table_file == '' or output_table_file == '':
        usage()
        sys.exit()
    
    # Read data table
    tb = Table.read(input_table_file, format='ascii')
    
    # Check data table row number
    if len(tb) >= 2:
        multi_row = True
    else:
        multi_row = False
    
    # loop each row
    for irow in range(len(tb)):
        ID = -99
        z = -99
        wave = []
        flux = []
        flux_err = []
        filter_name = []
        for icol, colname in enumerate(tb.columns):
            if icol == 0:
                ID = tb[colname][irow]
            elif icol == 1:
                if np.isnan(input_redshift):
                    z = tb[colname][irow]
                else:
                    z = input_redshift
            elif icol%2 == 0:
                flux.append(tb[colname][irow])
                filter_name.append(colname)
                wave.append(-99)
            elif icol%2 == 1:
                flux_err.append(tb[colname][irow])
        # 
        data_output = {}
        data_output['wave'] = wave
        data_output['flux'] = flux
        data_output['flux_err'] = flux_err
        data_output['unit'] = ['mJy']*len(flux)
        data_output['filter_name'] = filter_name
        tbout = Table(data_output)
        # 
        if multi_row:
            output_prefix = 'ID_%d/'%(ID)
            os.makedirs(output_prefix)
        else:
            output_prefix = ''
        tbout.write(output_prefix+output_table_file, format='ascii.fixed_width', delimiter='  ', bookend=True, overwrite=True)
        with open(output_prefix+output_table_file, 'r+') as fp:
            fp.seek(0)
            fp.write('#')
        print('Output to "%s"!'%(output_prefix+output_table_file))















