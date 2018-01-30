#!/usr/bin/env python
# 
# Aim: 
#      This code reads one solution of the michi2 SED fitting from the chisq table file, 
#      then goes to each SED library file which is used for the fitting, 
#      then extract all SED templates used for that fitting solution and also 
#      add all SED components together to make a total SED, 
#      and output to the specific directory given by the user. 
# Inputs:
#      First input is the michi2 SED fitting output chisq table file name. 
#      The second input should be a line number corresponding to the data rows (non-commented lines) in the first input chisq table. 
#      The third input is the output directory name. 
# 
# 

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')

from CrabTable import *

import numpy
import astropy
import astropy.io.ascii as asciitable
import scipy
import scipy.interpolate
from copy import copy




####################################
#               FUNC               #
####################################

def lib_file_get_header(Lib_file):
    # return the 'NVAR1' value and 
    # return the number of '#' commented lines at the beginning of the Lib_file.
    Lib_header = {}
    NVAR1 = 0
    NLINE = 0
    with open(Lib_file,'r') as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith('#'):
                NLINE = NLINE + 1
                if data_line.startswith('# NVAR1'):
                    data_line_split = data_line.split('=')
                    if len(data_line_split) >= 2:
                        data_line_split = data_line_split[1]
                        if data_line_split.find('#') > 0:
                            data_line_split = data_line_split.split('#')[0]
                        NVAR1 = int(data_line_split)
            else:
                break
        fp.close()
    Lib_header['NLINE'] = NLINE
    Lib_header['NVAR1'] = NVAR1
    return Lib_header


def lib_file_get_data_block(Lib_file, starting_data_line_index):
    # starting_data_line_index starts from 0 in the pure data block, i.e., no commented lines accounted. 
    # e.g., starting_data_line_index = 0, means the first SED template in the SED library file. 
    Lib_header = lib_file_get_header(Lib_file)
    Lib_begin = Lib_header['NLINE'] + starting_data_line_index                         # the line number index (starting from 0) in the Lib_file, which defines the data block of one SED template.
    Lib_end   = Lib_header['NLINE'] + starting_data_line_index + Lib_header['NVAR1']-1 # the line number index (starting from 0) in the Lib_file, which defines the data block of one SED template.
    print('numpy.genfromtxt(Lib_file, skip_header=%d, max_rows=%d)'%(Lib_begin, Lib_header['NVAR1']))
    Lib_arr = numpy.genfromtxt(Lib_file, skip_header=Lib_begin, max_rows=Lib_header['NVAR1'])
    return Lib_arr


def check_array(arr):
    # check array type numpy.ndarray
    if type(arr) is not list and type(arr) is not numpy.ndarray:
        arr = [arr]
    if type(arr) is list:
        arr = numpy.array(arr)
    return arr

def nearest_interp(xi, x, y):
    # https://stackoverflow.com/questions/21002799/extraploation-with-nearest-method-in-python
    idx = numpy.abs(x - xi[:,None])
    return y[idx.argmin(axis=1)]

def spline(input_x, input_y, output_x, xlog=0, ylog=0, outputxlog=None, outputylog=None, fill=numpy.nan):
    # spline
    # note that we can input xlog, ylog, outputxlog, outputylog to deal with logarithm input/output needs
    # xlog>0 means the input_x array is in linear space and will be converted to log space when doing the spline
    # ylog>0 means the input_y array is in linear space and will be converted to log space when doing the spline
    # outputxlog>0 means the output_x array is in linear space and will be converted to log space when doing the spline
    # outputylog>0 means the output_y array is in linear space but is in log space when doing the spline, so we need to convert it back to linear space after the spline.
    input_x_coord = check_array(input_x)
    input_y_value = check_array(input_y)
    output_x_coord = check_array(output_x)
    if xlog>0:
        input_x_mask = (input_x_coord<=0.0)
        input_x_coord[input_x_mask] = numpy.nan
        input_x_coord = numpy.log10(input_x)
    if ylog>0:
        input_y_mask = (input_y_value<=0.0)
        input_y_value[input_y_mask] = numpy.nan
        input_y_value = numpy.log10(input_y_value)
    # 
    if outputxlog is None:
        outputxlog = xlog
    if outputylog is None:
        outputylog = ylog
    # 
    if outputxlog>0:
        output_x_mask = (output_x_coord<=0.0)
        output_x_coord[output_x_mask] = numpy.nan
        output_x_coord = numpy.log10(output_x_coord)
    # 
    #print(numpy.column_stack((input_x_coord, input_y_value)))
    input_mask = numpy.isnan(input_y_value)
    
    #########output_y_value = scipy.interpolate.UnivariateSpline(input_x_coord, input_y_value, w=~input_mask)(output_x_coord)
    
    ###############output_y_value = scipy.interpolate.interpn(input_x_coord, input_y_value, output_x_coord, method='linear')
    
    ###########output_y_value = scipy.interpolate.interp1d(input_x_coord, input_y_value, kind='nearest')(output_x_coord) # error on outside data
    
    ###########output_y_value = numpy.interp() # this does the extrapolation, but see https://stackoverflow.com/questions/21002799/extraploation-with-nearest-method-in-python
    
    #output_y_value = nearest_interp(output_x_coord, input_x_coord, input_y_value)
    
    output_y_value = scipy.interpolate.spline(input_x_coord, input_y_value, output_x_coord, order='1') # order=3, kind='smoothest', conds=None
    
    #spl = scipy.interpolate.CubicSpline(input_x_coord, input_y_value) # axis=0, bc_type='not-a-knot', extrapolate=None
    #output_y_value = spl(output_x_coord)
    
    # 
    # deal with data out of X range
    output_mask = (output_x_coord<numpy.nanmin(input_x_coord)) | (output_x_coord>numpy.nanmax(input_x_coord))
    output_y_value[output_mask] = fill
    # 
    if outputylog>0:
        output_y = numpy.power(10,output_y_value)
        output_y[output_mask] = fill #<TODO># fill with 0.0 instead of nan
    else:
        output_y = output_y_value
    # 
    return output_y














####################################
#               MAIN               #
####################################

if not len(sys.argv) > 3:
    print('Usage: michi2_read_lib_SEDs.py chisq_file line_number output_directory')
    print('Example: michi2_read_lib_SEDs.py fit_5.out 1 output_SEDs_i_300')
    sys.exit()

chisq_file = sys.argv[1]
info_file = sys.argv[1]+'.info'
line_number = int(sys.argv[2])
output_dir = sys.argv[3]

if not os.path.isfile(chisq_file):
    print('Error! "%s" was not found!'%(chisq_file))
    sys.exit()

if not os.path.isfile(info_file):
    print('Error! "%s" was not found!'%(info_file))
    sys.exit()

if not (line_number>0):
    print('Error! The input line numebr %s is non-positive!')
    sys.exit()

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)



# Read chisq table
chisq_table = CrabTable(chisq_file, verbose=0)

# Fix data table header problem
chisq_table_headers = []
with open(chisq_file,'r') as fp:
    while True:
        data_line = fp.readline()
        if not data_line:
            break
        if data_line.startswith('#'):
            data_line_split = data_line.replace('#','').strip().split()
            if len(data_line_split) == len(chisq_table.TableHeaders):
                chisq_table_headers = data_line_split
        else:
            break
    fp.close()

# Read info table
info_table = CrabTableReadInfo(info_file, verbose=0)

#print(chisq_table.TableHeaders)
#print(chisq_table_headers)

# Loop each SED library which are listed in the info table and have columns in the chisq table.
Lib_number = int(info_table['NLIB'])
Lib_array = {}
Lib_array['TOT'] = {}
Lib_array['TOT']['log_X'] = numpy.arange(-2,6,0.001) # wavelength_um grid
Lib_array['TOT']['X'] = numpy.power(10, Lib_array['TOT']['log_X']) # make it in linear space
Lib_array['TOT']['Y'] = Lib_array['TOT']['X'] * 0.0
for iLib in range(Lib_number):
    # check SED lib file
    Lib_file = info_table['LIB%d'%(iLib+1)]
    if not os.path.isfile(Lib_file):
        print('Error! "%s" was not found!'%(Lib_file))
        sys.exit()
    # 
    # get one SED lib template according to the user input 'line_number'
    Lib_icol = 2+2*iLib
    Lib_acol = 2+2*iLib+1
    Lib_istr = chisq_table_headers[Lib_icol]
    Lib_astr = chisq_table_headers[Lib_acol]
    Lib_i = chisq_table.getColumn(Lib_icol+1)
    Lib_a = chisq_table.getColumn(Lib_acol+1)
    print(Lib_istr, Lib_i[line_number-1])
    print(Lib_astr, Lib_a[line_number-1])
    #print(chisq_table.getColumn(Lib_icol))
    os.system('echo "%s" > "%s/line_number"'%(line_number, output_dir))
    os.system('echo "%s" > "%s/%s"'%(Lib_i[line_number-1], output_dir, Lib_istr))
    os.system('echo "%s" > "%s/%s"'%(Lib_a[line_number-1], output_dir, Lib_astr))
    # 
    # read lib data block from line file, starting from the data line index 'Lib_i[line_number-1]'
    Lib_arr = lib_file_get_data_block(Lib_file, Lib_i[line_number-1])
    Lib_x = Lib_arr[:,0]
    Lib_y = Lib_arr[:,1] * Lib_a[line_number-1]
    Out_file = output_dir+os.sep+'SED_LIB%d'%(iLib+1)
    asciitable.write(numpy.column_stack((Lib_x,Lib_y)), 
                        Out_file, 
                        Writer=asciitable.FixedWidthTwoLine, 
                        names=['X', 'Y'], 
                        overwrite=True)
    os.system('sed -i.bak -e "1s/^/#/" "%s"'%(Out_file))
    os.system('sed -i.bak -e "2s/^[ -]*/#/" "%s"'%(Out_file))
    print('Output to "%s"'%(Out_file))
    # 
    # sum to make total SED (only when a>0.0)
    if Lib_a[line_number-1] > 0.0:
        Lib_array['LIB%d'%(iLib+1)] = {}
        Lib_array['LIB%d'%(iLib+1)]['X'] = Lib_x
        Lib_array['LIB%d'%(iLib+1)]['Y'] = Lib_y
        Lib_array['LIB%d'%(iLib+1)]['log_X'] = numpy.log10(Lib_x)
        Lib_array['LIB%d'%(iLib+1)]['log_Y'] = numpy.log10(Lib_y)
        #print(Lib_array['LIB%d'%(iLib+1)]['log_X'])
        #print(Lib_array['LIB%d'%(iLib+1)]['log_Y'])
        Lib_array['TOT']['Y'] = Lib_array['TOT']['Y'] + spline(Lib_array['LIB%d'%(iLib+1)]['X'], Lib_array['LIB%d'%(iLib+1)]['Y'], Lib_array['TOT']['X'], xlog=1, ylog=1, fill=0.0)

Out_file = output_dir+os.sep+'SED_SUM'
asciitable.write(numpy.column_stack((Lib_array['TOT']['X'],Lib_array['TOT']['Y'])), 
                    Out_file, 
                    Writer=asciitable.FixedWidthTwoLine, 
                    names=['X', 'Y'], 
                    overwrite=True)
os.system('sed -i.bak -e "1s/^/#/" "%s"'%(Out_file))
os.system('sed -i.bak -e "2s/^[ -]*/#/" "%s"'%(Out_file))
print('Output to "%s"'%(output_dir+os.sep+'SED_SUM'))

#if not os.path.isfile('obj_%d/SED_LIB%d'%(i+1,j+1)):
#    BashCommand = 'cd obj_%d/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_SED ../%s %d %s SED_LIB%d'%\
#                        (i+1, \
#                            InfoDict['LIB%d'%(j+1)], \
#                                DataArray['i%d'%(j+1)][SelectIndex[i]], \
#                                    DataArray['a%d'%(j+1)][SelectIndex[i]], \
#                                        j+1)
#    print(BashCommand)
#    os.system(BashCommand)





