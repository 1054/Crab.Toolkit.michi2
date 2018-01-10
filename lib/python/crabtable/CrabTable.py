#!/usr/bin/env python2.7
# 

################################
# 
# class CrabTable
# 
#   Example: 
#            DataTable = CrabTable()
#            DataTable.load('aaa.fits')
#            DataTable.load('aaa.ascii')
# 
#   Last update: 
#            20170930 12h23m -- copied from '/Users/dzliu/Cloud/Github/Crab.Toolkit.CAAP/bin/CrabDataTable.py'
# 
################################

try:
    import pkg_resources
except ImportError:
    raise SystemExit("Error! Failed to import pkg_resources!")

pkg_resources.require("numpy")
pkg_resources.require("astropy>=1.3")

import os
import sys
import glob
import math
import numpy
import astropy
from astropy import units
from astropy.io import fits
import astropy.io.ascii as asciitable













# 
class CrabTable(object):
    # 
    def __init__(self, data_table='', fits_extension=0, fix_string_columns=0, verbose=1):
        self.DataTableFile = ''
        self.DataTableFormat = ''
        self.DataTableStruct = []
        self.DataTableIndex = []
        self.TableData = []
        self.TableHeaders = []
        self.TableColumns = []
        self.TableIndex = 0
        self.World = {}
        self.World['verbose'] = verbose
        if data_table != '':
            self.load(data_table, fits_extension, fix_string_columns)
        # 
        #print a column
        #print(self.TableData.field('FWHM_MAJ_FIT'))
    # 
    def clear(self):
        self.DataTableFile = ''
        self.DataTableFormat = ''
        self.DataTableStruct = []
        self.DataTableIndex = []
        self.TableData = []
        self.TableHeaders = []
        self.TableColumns = []
        self.TableIndex = 0
    # 
    def load(self, data_table, fits_extension=0, fix_string_columns=0): 
        self.clear()
        if data_table != '':
            # get file path
            self.DataTableFile = data_table
            # check file format
            if data_table.endswith('.fits') or data_table.endswith('.FITS'):
                self.DataTableFormat = 'FITS'
            else:
                self.DataTableFormat = 'ASCII'
            # print info
            if not 'verbose' in self.World:
                self.World['verbose'] = 1
            if self.World['verbose'] > 0:
                print('Reading %s Table: %s'%(self.DataTableFormat, self.DataTableFile))
            # read data table
            if self.DataTableFormat == 'FITS':
                # open FITS format data table with astropy.io.fits
                self.DataTableStruct = fits.open(self.DataTableFile)
                #print TableStruct.info()
                # get the fits extension id according to the input fits_extension
                # skip astropy.io.fits.hdu.image.PrimaryHDU
                # only consider astropy.io.fits.hdu.table.BinTableHDU
                self.DataTableIndex = []
                for i in range(len(self.DataTableStruct)):
                    #print(type(self.DataTableStruct[i]))
                    if type(self.DataTableStruct[i]) is astropy.io.fits.hdu.table.BinTableHDU:
                        self.DataTableIndex.append(i)
                # 
                if len(self.DataTableIndex)>0:
                    self.TableIndex = fits_extension
                    self.TableData = self.DataTableStruct[self.DataTableIndex[self.TableIndex]].data
                    self.TableColumns = self.DataTableStruct[self.DataTableIndex[self.TableIndex]].columns # dtype TableColumns
                    self.TableHeaders = self.TableColumns.names
                    # deal with data units
                    #DataTable.DataTableStruct[1].header['TUNIT455']
                    # deal with string column dtype, make sure string column have string width at least 100 chars.
                    if fix_string_columns > 0:
                        for i in range(len(self.TableHeaders)):
                            if self.TableData[self.TableHeaders[i]].dtype.char == 'S':
                                #print((self.TableData.dtype))
                                dt = self.TableData.dtype
                                ds = dt.descr
                                dw = int(ds[i][1][ds[i][1].index('S')+1::])
                                if dw < 100:
                                    ds[i] = (ds[i][0], 'S100')
                                    dt = numpy.dtype(ds)
                                    self.TableData = self.TableData.astype(dt)
                                #print((self.TableData.dtype))
                                # -- https://stackoverflow.com/questions/9949427/how-to-change-the-dtype-of-a-numpy-recarray
                                #TempArray = numpy.array(self.TableData[self.TableHeaders[i]], dtype=object)
                                #TempColumn = astropy.io.fits.Column(format='A100', array=TempArray)
                                #self.TableData[self.TableHeaders[i]] = TempColumn
                                ##self.TableData[self.TableHeaders[i]] = TempArray
                                ##print(type(self.TableData))
                                ##print(type(self.TableData.array)) # not working
                                ##print(self.TableData._coldefs)
                                #print(self.TableData._coldefs[self.TableHeaders[i]])
                                #print(type(self.TableData._coldefs[self.TableHeaders[i]]))
                                ##self.TableData._coldefs[self.TableHeaders[i]].format.dtype = object
                                ##print(self.TableData._coldefs[self.TableHeaders[i]])
                                ##print(type(self.TableData._coldefs[self.TableHeaders[i]]))
                                #print((self.TableData.dtype))
                                # -- http://docs.astropy.org/en/stable/_modules/astropy/io/fits/fitsrec.html
                                # -- http://docs.astropy.org/en/stable/_modules/astropy/io/fits/column.html
                else:
                    print('Error! Could not find astropy.io.fits.hdu.table.BinTableHDU from the fits file %s with fits extension %d!'%(self.DataTableFile, fits_extension))
                    self.World['status'] = 'failed to load data table'
                    return
            else:
                # open ASCII format data table with astropy.io.ascii
                # http://cxc.harvard.edu/contrib/asciitable/
                self.DataTableStruct = asciitable.read(self.DataTableFile)
                self.TableIndex = 0
                self.TableData = self.DataTableStruct
                self.TableColumns = self.DataTableStruct.columns # dtype TableColumns
                self.TableHeaders = self.DataTableStruct.colnames
                # deal with string column dtype
                if fix_string_columns > 0:
                    for i in range(len(self.TableHeaders)):
                        if self.TableData[self.TableHeaders[i]].dtype.char == 'S' or self.TableData[self.TableHeaders[i]].dtype is numpy.string_:
                            TempArray = numpy.array(self.TableData[self.TableHeaders[i]], dtype=object)
                            self.TableData[self.TableHeaders[i]] = TempArray
        else:
            print('Error! The input data table is an empty string!')
    # 
    def getData(self):
        return self.TableData
    # 
    def getRowNumber(self):
        return len(self.TableData)
    # 
    def getColumnNames(self):
        return self.TableHeaders
    # 
    def getColumn(self, ColNameOrNumb):
        if type(ColNameOrNumb) is str or type(ColNameOrNumb) is numpy.string_:
            if ColNameOrNumb in self.TableHeaders:
                return self.TableData.field(ColNameOrNumb)
            else:
                print("Error! Column name \"%s\" was not found in the data table!"%(ColNameOrNumb))
                return []
        else:
            if ColNameOrNumb >= 1 and ColNameOrNumb <= len(self.TableHeaders):
                return self.TableData.field(self.TableHeaders[int(ColNameOrNumb)-1])
            else:
                print("Error! Column number %d is out of allowed range (1 - %d)!"%(int(ColNameOrNumb),len(self.TableHeaders)))
                return []
    # 
    def getColumnIndex(self, ColNameOrNumb):
        if type(ColNameOrNumb) is str or type(ColNameOrNumb) is numpy.string_:
            if ColNameOrNumb in self.TableHeaders:
                for i in range(len(self.TableHeaders)):
                    if self.TableHeaders[i] == ColNameOrNumb:
                        return i
            else:
                print("Error! Column name \"%s\" was not found in the data table!"%(ColNameOrNumb))
                return []
        else:
            if ColNameOrNumb >= 1 and ColNameOrNumb <= len(self.TableHeaders):
                return ColNameOrNumb-1
            else:
                print("Error! Column number %d is out of allowed range (1 - %d)!"%(int(ColNameOrNumb),len(self.TableHeaders)))
                return []
    # 
    def setColumn(self, ColNameOrNumb, DataArray):
        if type(ColNameOrNumb) is str or type(ColNameOrNumb) is numpy.string_:
            if ColNameOrNumb in self.TableHeaders:
                self.TableData[ColNameOrNumb] = DataArray
            else:
                print("Error! Column name \"%s\" was not found in the data table!"%(ColNameOrNumb))
                return
        else:
            if ColNameOrNumb >= 1 and ColNameOrNumb <= len(self.TableHeaders):
                self.TableData[self.TableHeaders[int(ColNameOrNumb)-1]] = DataArray
            else:
                print("Error! Column number %d is out of allowed range (1 - %d)!"%(int(ColNameOrNumb),len(self.TableHeaders)))
                return
    # 
    def setCell(self, ColNameOrNumb, RowNumb, DataValue):
        # ColNumb starts from 1
        # RowNumb starts from 1
        if type(ColNameOrNumb) is str or type(ColNameOrNumb) is numpy.string_:
            if ColNameOrNumb in self.TableHeaders:
                ColName = ColNameOrNumb
                if RowNumb >= 1 and RowNumb <= len(self.TableData[ColName]):
                    self.TableData[ColName][int(RowNumb)-1] = DataValue
                else:
                    print("Error! Row number %d is out of allowed range (1 - %d)!"%(int(RowNumb),len(self.TableData[ColName])))
            else:
                print("Error! Column name \"%s\" was not found in the data table!"%(ColNameOrNumb))
                return
        else:
            if ColNameOrNumb >= 1 and ColNameOrNumb <= len(self.TableHeaders):
                ColName = self.TableHeaders[int(ColNameOrNumb)-1]
                if RowNumb >= 1 and RowNumb <= len(self.TableData[ColName]):
                    self.TableData[ColName][int(RowNumb)-1] = DataValue
                else:
                    print("Error! Row number %d is out of allowed range (1 - %d)!"%(int(RowNumb),len(self.TableData[ColName])))
            else:
                print("Error! Column number %d is out of allowed range (1 - %d)!"%(int(ColNameOrNumb),len(self.TableHeaders)))
                return
    # 
    def saveAs(self, OutputFilePath, overwrite=False, writer=asciitable.FixedWidthTwoLine):
        # backup output file
        if self.DataTableStruct:
            has_overwritten = False
            if os.path.isfile(OutputFilePath):
                if overwrite == True:
                    os.system("mv %s %s"%(OutputFilePath, OutputFilePath+'.backup'))
                    has_overwritten = True
                else:
                    print("We will not overwrite unless you specify saveAs(overwrite=True)!")
                    return
            # 
            if OutputFilePath.endswith('.fits') or OutputFilePath.endswith('.FITS'):
                #<TODO># self.DataTableStruct[self.DataTableIndex[self.TableIndex]].data = self.TableData
                self.DataTableStruct.writeto(OutputFilePath)
            else:
                #self.DataTableStruct = self.TableData
                asciitable.write(self.DataTableStruct, OutputFilePath, Writer=writer)
            # print output file path
            if has_overwritten:
                print("Output to %s! (A backup has been created as %s)"%(OutputFilePath, OutputFilePath+'.backup'))
            else:
                print("Output to %s!"%(OutputFilePath))











# 
def CrabTableReadInfo(data_table, fits_extension=0, key_name=[], verbose=1):
    # This function reads all rows in a text file or data table with content like
    # aaa = aaa_value
    # bbb = bbb_value
    # ccc = ccc_value
    # 
    # prepare output info dict
    InfoDict = {}
    # 
    # read input data table
    if os.path.isfile(data_table):
        if data_table.endswith('.fits') or data_table.endswith('.FITS'):
            data_table_format = 'FITS'
        else:
            data_table_format = 'ASCII'
        # print message
        if verbose > 0:
            print('Reading %s Table: %s'%(data_table_format, data_table))
        # read data table
        if data_table_format == 'FITS':
            # open FITS format data table with astropy.io.fits
            data_table_struct = fits.open(data_table)
            #print TableStruct.info()
            # get the fits extension id according to the input fits_extension
            # skip astropy.io.fits.hdu.image.PrimaryHDU
            # only consider astropy.io.fits.hdu.table.BinTableHDU
            data_table_index = []
            for i in range(len(data_table_struct)):
                if type(data_table_struct[i]) is astropy.io.fits.hdu.table.BinTableHDU:
                    data_table_index.append(i)
            # 
            if len(data_table_index)>fits_extension:
                TableIndex = data_table_index[fits_extension]
                TableData = data_table_struct[data_table_index[fits_extension]].data
                TableColumns = data_table_struct[data_table_index[fits_extension]].columns # dtype TableColumns
                TableHeaders = TableColumns.names
            # loop each table row, create dict data
            if len(TableHeaders) > 1:
                for i in range(len(TableData)):
                    InfoDict[TableData[TableHeaders[0]]][i] = TableData[TableHeaders[1]][i]
        else:
            # open ASCII format data table with astropy.io.ascii
            # http://cxc.harvard.edu/contrib/asciitable/
            data_table
            with open(data_table, "r") as data_table_ptr:
                data_table_lines = data_table_ptr.readlines()
                for data_table_line in data_table_lines:
                    data_table_line_items = data_table_line.split('=')
                    if len(data_table_line_items)>=2:
                        InfoDict[data_table_line_items[0].strip()] = data_table_line_items[1].strip()
    # 
    # return
    return InfoDict




























