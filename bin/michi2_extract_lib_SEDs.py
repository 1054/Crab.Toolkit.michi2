#!/usr/bin/env python
# 
# Aim: 
#      This code reads the michi2 SED LIB file and extract a single SED LIB acoording to the user-given index. 
# Inputs:
#      First input is the michi2 SED LIB file name. 
#      The second input should be the line index corresponding to the LIB data block rows (non-commented lines). 
#      The third argument is the output single LIB file name. 
# 
# 

import os
import sys
import re
import time

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

def lib_header_get_key(Lib_header, Key_name):
    Key_value = ''
    Key_comment = ''
    if type(Lib_header) is dict:
        if 'LINES' in Lib_header:
            if len(Lib_header['LINES']) > 0:
                for header_line in Lib_header['LINES']:
                    matcher = re.match('^# *%s *= *([^#]+)( *#* *)(.*)$'%(Key_name), header_line)
                    if matcher:
                        Key_value = matcher.group(1).strip()
                        print('Key_value = %s'%(Key_value))
                        if len(matcher.groups()) >= 3:
                            Key_comment = matcher.group(3)
                            print('Key_comment = %s'%(Key_comment))
    return Key_value, Key_comment

def lib_header_mod_key(Lib_header, Key_name, Key_value, Key_comment=None):
    Key_modified = False
    if type(Lib_header) is dict:
        if 'LINES' in Lib_header:
            if len(Lib_header['LINES']) > 0:
                for i in range(len(Lib_header['LINES'])):
                    header_line = Lib_header['LINES'][i]
                    matcher = re.match('^# *%s *= *([^#]+)( *#* *)(.*)$'%(Key_name), header_line)
                    if matcher:
                        # found the Key_name
                        if Key_comment is None:
                            # try to get existing Key_comment
                            if len(matcher.groups()) >= 3:
                                Key_comment = matcher.group(3)
                        # update Lib_header
                        if Key_comment is not None:
                            Lib_header['LINES'][i] = '# %s = %s # %s\n'%(Key_name, Key_value, Key_comment)
                        else:
                             Lib_header['LINES'][i] = '# %s = %s\n'%(Key_name, Key_value) 
                        Key_modified = True
    return Key_modified

def lib_file_get_header(Lib_file):
    # return the 'NVAR1' value and 
    # return the number of '#' commented lines at the beginning of the Lib_file.
    Lib_header = {}
    Lib_header['LINES'] = []
    NVAR1 = 0
    NLINE = 0
    with open(Lib_file,'r') as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith('#'):
                NLINE = NLINE + 1
                Lib_header['LINES'].append(data_line)
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


def lib_file_get_data_block(Lib_file, starting_data_line_index, Lib_header = []):
    # starting_data_line_index starts from 0 in the pure data block, i.e., no commented lines accounted. 
    # e.g., starting_data_line_index = 0, means the first SED template in the SED library file. 
    if Lib_header == []: Lib_header = lib_file_get_header(Lib_file)
    Lib_begin = Lib_header['NLINE'] + starting_data_line_index                         # the line number index (starting from 0) in the Lib_file, which defines the data block of one SED template.
    Lib_end   = Lib_header['NLINE'] + starting_data_line_index + Lib_header['NVAR1']-1 # the line number index (starting from 0) in the Lib_file, which defines the data block of one SED template.
    print('numpy.genfromtxt(Lib_file, skip_header=%d, max_rows=%d)'%(Lib_begin, Lib_header['NVAR1']))
    Lib_arr = numpy.genfromtxt(Lib_file, skip_header=Lib_begin, max_rows=Lib_header['NVAR1'])
    return Lib_arr


def lib_file_get_data_block_lines(Lib_file, starting_data_line_index, Lib_header = []):
    # starting_data_line_index starts from 0 in the pure data block, i.e., no commented lines accounted. 
    # e.g., starting_data_line_index = 0, means the first SED template in the SED library file. 
    if Lib_header == []: Lib_header = lib_file_get_header(Lib_file)
    Lib_begin = Lib_header['NLINE'] + starting_data_line_index                         # the line number index (starting from 0) in the Lib_file, which defines the data block of one SED template.
    Lib_end   = Lib_header['NLINE'] + starting_data_line_index + Lib_header['NVAR1']-1 # the line number index (starting from 0) in the Lib_file, which defines the data block of one SED template.
    Lib_data_block_lines = []
    with open(Lib_file) as fp:
        for i, line in enumerate(fp):
            if i < Lib_begin:
                pass
            elif i <= Lib_end:
                Lib_data_block_lines.append(line)
            elif i > Lib_end:
                break
    return Lib_data_block_lines














####################################
#               MAIN               #
####################################

if not len(sys.argv) > 3:
    print('Usage: michi2_extract_lib_SEDs.py LIB_file LIB_index output_fileectory')
    print('Example: michi2_extract_lib_SEDs.py lib.full.SED 0 output.lib.single.SED')
    sys.exit()

LIB_file = sys.argv[1]
LIB_index = int(sys.argv[2])
output_file = sys.argv[3]

if not os.path.isfile(LIB_file):
    print('Error! "%s" was not found!'%(LIB_file))
    sys.exit()

if not (LIB_index>=0):
    print('Error! The input line numebr %s is negative!')
    sys.exit()





# Read SED LIB file
LIB_header = lib_file_get_header(LIB_file)
LIB_data_block_lines = lib_file_get_data_block_lines(LIB_file, LIB_index, Lib_header=LIB_header)

NPAR = int(lib_header_get_key(LIB_header, 'NPAR')[0])

i=1
while i <= NPAR:
    lib_header_get_key(LIB_header, 'NPAR%d'%(i))
    i = i + 1
lib_header_get_key(LIB_header, 'NVAR2')

i=1
while i <= NPAR:
    lib_header_mod_key(LIB_header, 'NPAR%d'%(i), 1)
    i = i + 1
lib_header_mod_key(LIB_header, 'NVAR2', 1)

i=1
while i <= NPAR:
    lib_header_get_key(LIB_header, 'NPAR%d'%(i))
    i = i + 1
lib_header_get_key(LIB_header, 'NVAR2')




# Update current time
LIB_header['LINES'][0] = '# %s (extracted from \"%s\" at index %d)\n'%(time.strftime('%Y-%m-%d %H:%M:%S %Z'), LIB_file, LIB_index)



# Output extracted single LIB file
with open(output_file, 'w') as fp:
    for header_line in LIB_header['LINES']:
        fp.write(header_line)
    for i in range(len(LIB_data_block_lines)):
        fp.write(LIB_data_block_lines[i])
    print('Output to "%s"!'%(output_file))







