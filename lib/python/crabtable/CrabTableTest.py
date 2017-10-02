#!/usr/bin/env python2.7
# 

import os, sys
from CrabTable import *
from pprint import pprint

datatable_ascii = CrabTable(os.path.dirname(sys.argv[0])+os.sep+'CrabTableTest.ascii')
datatable_fits = CrabTable(os.path.dirname(sys.argv[0])+os.sep+'CrabTableTest.fits')


#asciitable.write(datatable_ascii.TableData, sys.stdout, Writer=asciitable.FixedWidthTwoLine)
#asciitable.write(datatable_fits.TableData, sys.stdout, Writer=asciitable.FixedWidthTwoLine)
pprint(datatable_ascii.TableData)
pprint(datatable_fits.TableData)


datatable_ascii.setCell(1,1,99.0)
datatable_fits.setCell(1,1,99.0)

datatable_ascii.setCell('str',3,'test')
datatable_fits.setCell('str',3,'test')


pprint(datatable_ascii.TableData)
pprint(datatable_fits.TableData)


datatable_ascii.saveAs(os.path.dirname(sys.argv[0])+os.sep+'CrabTableTest.ascii.2.ascii', overwrite=True)
datatable_fits.saveAs(os.path.dirname(sys.argv[0])+os.sep+'CrabTableTest.fits.2.fits', overwrite=True)





