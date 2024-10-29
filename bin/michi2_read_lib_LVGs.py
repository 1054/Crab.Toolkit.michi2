#!/usr/bin/env python
# 
# Aim: 
#      This code reads one solution of the michi2 LVG fitting from the chisq table file, 
#      then goes to each LVG library file which is used for the fitting, 
#      then extract all LVG templates used for that fitting solution and also 
#      add all LVG components together to make a total LVG, 
#      and output to the specific directory given by the user. 
# Inputs:
#      First input is the michi2 LVG fitting output chisq table file name. 
#      The second input should be a line number corresponding to the data rows (non-commented lines) in the first input chisq table. 
#      The third input is the output directory name. 
# 
# 

import os, sys, re

sys.path.append(os.path.dirname(os.path.dirname(__file__)) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')

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
    NLINE = 0
    NBYTE = 0
    with open(Lib_file,'r') as fp:
        last_valid_header_line = ''
        first_valid_data_line = ''
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith('#'):
                NLINE = NLINE + 1
                #NBYTE = fp.tell() # the last +1 accounts for the ending 
                if data_line.strip() != '#':
                    for header_key in ['NVAR1', 'NVAR2', \
                                       'NPAR', \
                                       'NPAR1', 'NPAR2', 'NPAR3', 'NPAR4', 'NPAR5', 'NPAR6', 'NPAR7', 'NPAR8', 'NPAR9', 'NPAR10', 'NPAR11', 'NPAR12', 'NPAR13', 'NPAR14', 'NPAR15', 'NPAR16', 'NPAR17', 'NPAR18', 'NPAR19', 'NPAR20', \
                                       'CPAR1', 'CPAR2', 'CPAR3', 'CPAR4', 'CPAR5', 'CPAR6', 'CPAR7', 'CPAR8', 'CPAR9', 'CPAR10', 'CPAR11', 'CPAR12', 'CPAR13', 'CPAR14', 'CPAR15', 'CPAR16', 'CPAR17', 'CPAR18', 'CPAR19', 'CPAR20', \
                                       'TPAR1', 'TPAR2', 'TPAR3', 'TPAR4', 'TPAR5', 'TPAR6', 'TPAR7', 'TPAR8', 'TPAR9', 'TPAR10', 'TPAR11', 'TPAR12', 'TPAR13', 'TPAR14', 'TPAR15', 'TPAR16', 'TPAR17', 'TPAR18', 'TPAR19', 'TPAR20', ]:
                        #<TODO># max 20 pars allowed.
                        #if data_line.startswith('# %s'%(header_key)):
                        if re.match(r'^#\s*%s\s*=.*'%(header_key), data_line):
                            data_line_split = data_line.split('=')
                            if len(data_line_split) >= 2:
                                data_line_split = data_line_split[1]
                                if data_line_split.find('#') > 0:
                                    data_line_split = data_line_split.split('#')[0]
                                if header_key.startswith('N') or header_key.startswith('C'):
                                    Lib_header[header_key] = int(data_line_split)
                                else:
                                    Lib_header[header_key] = data_line_split.strip()
                    last_valid_header_line = data_line
            elif data_line.strip() == '':
                continue
            else:
                first_valid_data_line = data_line
                NBYTE = fp.tell() - len(data_line)
                break
        fp.close()
    Lib_header['NLINE'] = NLINE
    Lib_header['NBYTE'] = NBYTE
    # 
    if last_valid_header_line == '':
        raise Exception('Error! Failed to read header column names from "%s"!'%(Lib_file))
    else:
        Lib_header['colnames'] = last_valid_header_line.replace('#','').strip().split()
    # 
    if len(Lib_header['colnames']) != Lib_header['NPAR'] + 2:
        raise Exception('Error! The column number and header NPAR do not match in "%s"!'%(Lib_file))
    # 
    if first_valid_data_line == '':
        raise Exception('Error! Failed to read a first data line from "%s"!'%(Lib_file))
    else:
        data_array = first_valid_data_line.strip().split()
        Lib_header['coltypes'] = []
        for i,data_str in enumerate(data_array):
            if re.match(r'^[0-9eE.+-]+$', data_str):
                Lib_header['coltypes'].append(numpy.float64)
                #Lib_header['coltypes'].append( (Lib_header['colnames'][i], float) )
            else:
                Lib_header['coltypes'].append(str) # numpy.str is a deprecated alias for the builtin `str`.
                #Lib_header['coltypes'].append( (Lib_header['colnames'][i], str) )
    # 
    return Lib_header


def lib_file_get_data_block(Lib_file, starting_data_line_index, Lib_header = [], verbose = True):
    # starting_data_line_index starts from 0 in the pure data block, i.e., no commented lines accounted. 
    # e.g., starting_data_line_index = 0, means the first LVG template in the LVG library file. 
    if Lib_header == []: Lib_header = lib_file_get_header(Lib_file)
    Lib_begin = Lib_header['NLINE'] + starting_data_line_index                         # the line number index (starting from 0) in the Lib_file, which defines the data block of one LVG template.
    Lib_end   = Lib_header['NLINE'] + starting_data_line_index + Lib_header['NVAR1']-1 # the line number index (starting from 0) in the Lib_file, which defines the data block of one LVG template.
    if verbose:
        print('Lib_begin %s'%(Lib_begin+1))
        print('Lib_end %s'%(Lib_end+1))
    # 
    #if verbose: print('numpy.genfromtxt(Lib_file, skip_header=%d, max_rows=%d)'%(Lib_begin, Lib_header['NVAR1']))
    #Lib_arr = numpy.genfromtxt(Lib_file, skip_header=Lib_begin, max_rows=Lib_header['NVAR1']) #<20200416># text will be lost as default dtype is float
    # 
    Lib_content = []
    with open(Lib_file, 'r') as fp:
        lines_to_skip = Lib_begin
        while lines_to_skip > 0:
            fp.readline()
            lines_to_skip -= 1
        lines_to_read = Lib_header['NVAR1']
        while lines_to_read > 0:
            Lib_content.append(tuple(fp.readline().strip().split()))
            lines_to_read -= 1
    # 
    Lib_arr = numpy.array(Lib_content)
    arr_dtypes = [('X','f8'), ('Y','f8')]
    for k in range(Lib_header['NPAR']):
        try:
            Lib_arr[0,k+2].astype(float)
            arr_dtypes.append(('PAR%d'%(k+1), 'f8'))
        except:
            max_len = numpy.max([len(t) for t in Lib_arr[:,k+2]])
            #arr_dtypes.append(('PAR%d'%(k+1), 'U%d'%(max_len+1)))
            arr_dtypes.append(('PAR%d'%(k+1), 'object'))
    # 
    Lib_arr = numpy.array(Lib_content, dtype=arr_dtypes)
    #print(Lib_arr.shape)
    #print(Lib_arr.dtype)
    # 
    return Lib_arr


def lib_file_get_data_block_quick(Lib_file, starting_data_line_index, every_data_line = 1, max_data_line_count = 0, Lib_header = [], verbose = True):
    # starting_data_line_index starts from 0 in the pure data block, i.e., no commented lines accounted. 
    # e.g., starting_data_line_index = 0, means the first LVG template in the LVG library file. 
    # 
    # https://stackoverflow.com/questions/19189961/python-fastest-access-to-line-in-file
    #    from itertools import islice
    #    from linecache import getline
    #    with open(Lib_file,'r') as fp:
    #        lines = list(islice(fp, 4003, 4005))
    #    for i in range(max_data_line_count)
    # 
    from itertools import islice
    if Lib_header == []: Lib_header = lib_file_get_header(Lib_file)
    if max_data_line_count == 0:
        #max_data_line_count = int(Lib_header['NVAR1']) * every_data_line #<BUG><20190920>#
        #max_data_line_count = int(Lib_header['NVAR2'])
        max_data_line_count = int(Lib_header['NVAR1'])
    Lib_arr = []
    with open(Lib_file,'r') as fp:
        # read lib file by slice
        Lib_lines = list(islice(fp, 
                                Lib_header['NLINE'] + starting_data_line_index, 
                                Lib_header['NLINE'] + starting_data_line_index + max_data_line_count * every_data_line, 
                                every_data_line))
        ##print('len(Lib_lines)', len(Lib_lines))
        ##print('Lib_lines[0]', Lib_lines[0], type(Lib_lines[0]))
        ##print('Lib_lines', ' '.join(Lib_lines))
        ##Lib_arr = numpy.fromstring(' '.join(Lib_lines), dtype=float, sep=' ')
        ##Lib_arr = numpy.fromstring(' '.join(Lib_lines), sep=' ') # dtype=float is removed because now we allow string value for some params
        #Lib_multiline_str = ' '.join([t.strip() for t in Lib_lines])
        #print('Lib_multiline_str', Lib_multiline_str)
        #Lib_arr = numpy.fromstring(Lib_multiline_str, sep=' ')
        #print('Lib_arr.shape', Lib_arr.shape)
        #print('Lib_arr', Lib_arr)
        ##print('reshape', len(Lib_lines), len(Lib_arr)/len(Lib_lines))
        ##Lib_arr = Lib_arr.reshape((len(Lib_lines), int(len(Lib_arr)/len(Lib_lines))) ) #<BUG><20190920>#
        #Lib_arr = Lib_arr.reshape((len(Lib_lines), int(len(Lib_arr)/len(Lib_lines))) )
        Lib_headline = '# X Y ' + ' '.join([Lib_header['TPAR%d'%(t)] for t in range(1,Lib_header['NPAR']+1)]) + '\n'
        Lib_table = CrabTable(data_lines=[Lib_headline]+Lib_lines, format='commented_header')
        Lib_arr = Lib_table.TableData
    # 
    #for iLib_line in range(max_data_line_count):
    #    Lib_arr.append( \
    #        numpy.fromstring( \
    #            linecache.getline(Lib_file, Lib_header['NLINE'] + starting_data_line_index + every_data_line * iLib_line), 
    #            dtype=float, sep=' '
    #        )
    #    ) # too slow
    #Lib_arr = numpy.array(Lib_arr)
    # 
    return Lib_arr


def lib_file_get_data_block_quick_2(Lib_file, starting_data_line_index, Lib_header, verbose = True):
    # starting_data_line_index starts from 0 in the pure data block, i.e., no commented lines accounted. 
    # e.g., starting_data_line_index = 0, means the first LVG template in the LVG library file. 
    import subprocess
    subproc = subprocess.run(['%s/michi2_read_lib_LVG'%(os.path.dirname(os.path.realpath(__file__))), Lib_file, '%d'%starting_data_line_index ], 
                             stdout=subprocess.PIPE)
    Flat_arr = numpy.fromstring(subproc.stdout, dtype=float, sep=' ')
    #print('len(Flat_arr) = %d'%(len(Flat_arr)))
    #print('size(Lib_arr) = (%s,%s)'%(Lib_header['NVAR1'], (len(Flat_arr)/Lib_header['NVAR1'])) )
    Lib_arr = Flat_arr.reshape( (int(Lib_header['NVAR1']), int(len(Flat_arr)/Lib_header['NVAR1'])) ) # numpy.column_stack((Flat_arr[::2],Flat_arr[1::2]))
    #print(Lib_arr[0])
    return Lib_arr


def lib_file_get_data_lines(Lib_file, select_nth_line, max_data_block_count = 0, Lib_header = [], verbose = True):
    # 20190920 new
    # select_nth_lines means the nth line in each data block (single model)
    # select_nth_lines can be a list or a single number
    # select_nth_lines 1 means the first line in each data block (LVG model)
    # we use mmap now
    # we allow multiple 'select_nth_lines' input as a list
    # 
    import shlex
    import mmap
    from itertools import islice
    # 
    if numpy.isscalar(select_nth_line):
        select_nth_lines = [select_nth_line]
    else:
        select_nth_lines = select_nth_line
    # 
    if Lib_header == []: 
        Lib_header = lib_file_get_header(Lib_file)
    # 
    model_rows = Lib_header['NVAR1'] # number of rows for each single model
    model_nums = Lib_header['NVAR2'] # number of models in the library
    model_cols = Lib_header['NPAR'] + 2 # number of parameters (columns) plus X and Y
    if max_data_block_count > 0:
        model_nums = max_data_block_count
    # 
    all_data_lines = numpy.full((len(select_nth_lines), model_nums, model_cols), '', dtype=object)
    #print('all_data_lines.shape', all_data_lines.shape, model_cols)
    # 
    with open(Lib_file, 'r') as fp:
        # memory-map the file, size 0 means whole file
        # note that for Windows OS use access=mmap.ACCESS_READ
        with mmap.mmap(fp.fileno(), 0, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ, offset=0) as mm:
            # skip header bytes
            mm.seek(Lib_header['NBYTE'])
            # read each model block
            data_lines = ['']*model_rows
            for i in range(model_rows * model_nums):
                data_line = mm.readline()
                #<debug>
                #if i == 0:
                #    print(data_line)
                for j in range(len(select_nth_lines)):
                    if select_nth_lines[j] == int(i % model_rows)+1:
                        if type(data_line) is bytes:
                            data_line = data_line.decode('utf-8')
                        #<debug>
                        #print(i, j, select_nth_lines[j], int(i % model_rows), int(i/model_rows))
                        #if data_line.find('C_atom')>=0:
                        #    print(i, j, select_nth_lines[j], int(i % model_rows), int(i/model_rows))
                        #    sys.exit()
                        all_data_lines[j,int(i/model_rows),:] = shlex.split(data_line)
    # 
    # convert str to float when needed
    #for k in range(all_data_lines.shape[2]):
    #    #print(k, all_data_lines[0,0,k])
    #    #if re.match(r'^[0-9]+$', all_data_lines[0,0,k]):
    #    #    all_data_lines[:,:,k] = all_data_lines[:,:,k].astype(numpy.int64)
    #    if re.match(r'^[0-9eE.+-]+$', all_data_lines[0,0,k]):
    #        all_data_lines[:,:,k] = all_data_lines[:,:,k].astype(numpy.float64)
    #    else:
    #        all_data_lines[:,:,k] = all_data_lines[:,:,k].astype(numpy.str)
    # 
    if numpy.isscalar(select_nth_line):
        return all_data_lines[0]
    else:
        return all_data_lines


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
    output_mask = numpy.logical_or(\
                        numpy.logical_or(\
                            output_x_coord<numpy.nanmin(input_x_coord), 
                            output_x_coord>numpy.nanmax(input_x_coord)
                        ), 
                        numpy.logical_or(\
                            numpy.isnan(output_y_value),
                            ~numpy.isfinite(output_y_value),
                        )
                    )
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

if __name__ == "__main__":
    
    if not len(sys.argv) > 3:
        print('Usage: michi2_read_lib_LVGs.py chisq_file line_number output_directory')
        print('Example: michi2_read_lib_LVGs.py fit_5.out 1 output_LVGs_i_300')
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
    
    # Read chisq table header from the first #-commented line (instead of the last one as asciitable does)
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
    
    # Loop each LVG library which are listed in the info table and have columns in the chisq table.
    Lib_number = int(info_table['NLIB'])
    Lib_array = {}
    Lib_array['TOT'] = {}
    Lib_array['TOT']['X'] = []
    Lib_array['TOT']['Y'] = []
    for iLib in range(Lib_number):
        # check LVG lib file
        Lib_file = info_table['LIB%d'%(iLib+1)]
        if not os.path.isfile(Lib_file):
            print('Error! "%s" was not found!'%(Lib_file))
            sys.exit()
        # 
        # get one LVG lib template according to the user input 'line_number'
        Lib_icol = 2+2*iLib
        Lib_acol = 2+2*iLib+1
        Lib_istr = chisq_table_headers[Lib_icol]
        Lib_astr = chisq_table_headers[Lib_acol]
        Lib_i = chisq_table.getColumn(Lib_icol+1)
        Lib_a = chisq_table.getColumn(Lib_acol+1)
        print(Lib_istr, Lib_i[line_number-1])
        print(Lib_astr, Lib_a[line_number-1])
        #print(chisq_table.getColumn(Lib_icol))
        Lib_header = lib_file_get_header(Lib_file)
        os.system('echo "%s" > "%s/line_number"'%(line_number, output_dir))
        os.system('echo "%s" > "%s/%s"'%(Lib_i[line_number-1], output_dir, Lib_istr))
        os.system('echo "%s" > "%s/%s"'%(Lib_a[line_number-1], output_dir, Lib_astr))
        os.system('echo "%d" > "%s/INDEX_LIB%d"'%(Lib_i[line_number-1]/int(Lib_header['NVAR1']), output_dir, iLib+1))
        # 
        # read lib data block from line file, starting from the data line index 'Lib_i[line_number-1]'
        Lib_arr = lib_file_get_data_block(Lib_file, Lib_i[line_number-1], Lib_header=Lib_header)
        Lib_x = Lib_arr['X']
        Lib_y = Lib_arr['Y'] * Lib_a[line_number-1]
        Out_file = output_dir+os.sep+'LVG_LIB%d'%(iLib+1)
        asciitable.write(numpy.column_stack((Lib_x,Lib_y)), 
                         Out_file, 
                         Writer=asciitable.FixedWidthTwoLine, 
                         names=['X', 'Y'], 
                         overwrite=True)
        os.system('sed -i.bak -e "1s/^/#/" "%s"'%(Out_file))
        os.system('sed -i.bak -e "2s/^[ -]*/#/" "%s"'%(Out_file))
        print('Output to "%s"'%(Out_file))
        # 
        # dump all lib columns
        Out_file = output_dir+os.sep+'LVG_LIB%d.par'%(iLib+1)
        Out_colnames = ['X', 'Y']
        for k in range(Lib_header['NPAR']):
            Out_colnames.append(Lib_header['TPAR%d'%(k+1)].replace('{','').replace('}','').replace('/','_'))
        asciitable.write(Lib_arr, 
                         Out_file, 
                         Writer=asciitable.FixedWidthTwoLine, 
                         names=Out_colnames, 
                         delimiter='  ', 
                         overwrite=True)
        os.system('sed -i.bak -e "1s/^/#/" "%s"'%(Out_file))
        os.system('sed -i.bak -e "2s/^[ -]*/#/" "%s"'%(Out_file))
        print('Output to "%s"'%(Out_file))
        # 
        # sum to make total LVG (only when a>0.0)
        if Lib_a[line_number-1] > 0.0:
            Lib_array['LIB%d'%(iLib+1)] = {}
            Lib_array['LIB%d'%(iLib+1)]['X'] = Lib_x
            Lib_array['LIB%d'%(iLib+1)]['Y'] = Lib_y
            Lib_array['LIB%d'%(iLib+1)]['log_X'] = numpy.log10(Lib_x)
            Lib_array['LIB%d'%(iLib+1)]['log_Y'] = numpy.log10(Lib_y)
            #print(Lib_array['LIB%d'%(iLib+1)]['log_X'])
            #print(Lib_array['LIB%d'%(iLib+1)]['log_Y'])
            if len(Lib_array['TOT']['X']) == 0:
                Lib_array['TOT']['X'] = Lib_x
                Lib_array['TOT']['Y'] = Lib_y
            else:
                Lib_array['TOT']['Y'] += Lib_y
    
    Out_file = output_dir+os.sep+'LVG_SUM'
    asciitable.write(numpy.column_stack((Lib_array['TOT']['X'],Lib_array['TOT']['Y'])), 
                        Out_file, 
                        Writer=asciitable.FixedWidthTwoLine, 
                        names=['X', 'Y'], 
                        overwrite=True)
    os.system('sed -i.bak -e "1s/^/#/" "%s"'%(Out_file))
    os.system('sed -i.bak -e "2s/^[ -]*/#/" "%s"'%(Out_file))
    print('Output to "%s"'%(output_dir+os.sep+'LVG_SUM'))
    
    #if not os.path.isfile('obj_%d/LVG_LIB%d'%(i+1,j+1)):
    #    BashCommand = 'cd obj_%d/; /Users/dzliu/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_read_lib_LVG ../%s %d %s LVG_LIB%d'%\
    #                        (i+1, \
    #                            InfoDict['LIB%d'%(j+1)], \
    #                                DataArray['i%d'%(j+1)][SelectIndex[i]], \
    #                                    DataArray['a%d'%(j+1)][SelectIndex[i]], \
    #                                        j+1)
    #    print(BashCommand)
    #    os.system(BashCommand)





