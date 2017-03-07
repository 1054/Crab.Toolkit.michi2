#!/usr/bin/env python2.7
# 

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
import re
import json



# 
# class CrabFitsTable
# 
# copied from $HOME/Cloud/Github/Crab.Toolkit/python/a_dzliu_python_lib_highz.py
# 
class CrabFitsTable(object):
    # 
    def __init__(self, FitsTableFile, FitsTableNumb=0):
        self.FitsTableFile = FitsTableFile
        print "# Reading Fits Table: %s"%(self.FitsTableFile)
        self.FitsStruct = fits.open(self.FitsTableFile)
        self.TableColumns = []
        self.TableData = []
        self.TableHeaders = []
        self.World = {}
        #print TableStruct.info()
        TableCount = 0
        for TableId in range(len(self.FitsStruct)):
            if type(self.FitsStruct[TableId]) is astropy.io.fits.hdu.table.BinTableHDU:
                if TableCount == FitsTableNumb:
                    self.TableColumns = self.FitsStruct[TableId].columns
                    self.TableData = self.FitsStruct[TableId].data
                TableCount = TableCount + 1
        if(TableCount==0):
            print "Error! The input FitsTableFile does not contain any data table!"
        else:
            self.TableHeaders = self.TableColumns.names
        #print a column
        #print(self.TableData.field('FWHM_MAJ_FIT'))
    # 
    def getData(self):
        return self.TableData
    # 
    def getColumnNames(self):
        return self.TableHeaders
    # 
    def getColumn(self, ColNameOrNumb):
        if type(ColNameOrNumb) is str:
            if ColNameOrNumb in self.TableHeaders:
                return self.TableData.field(ColNameOrNumb)
            else:
                print("Error! Column name \"%s\" was not found in the data table!"%(ColNameOrNumb))
                return []
        else:
            if ColNameOrNumb >= 0 and ColNameOrNumb < len(self.TableHeaders):
                return numpy.array(self.TableData[int(ColNameOrNumb)])
            else:
                print("Error! Column number %d is out of allowed range (0 - %d)!"%(int(ColNameOrNumb),len(self.TableHeaders)-1))
                return []
    # 
    def setColumn(self, ColNameOrNumb, DataArray):
        if type(ColNameOrNumb) is str:
            if ColNameOrNumb in self.TableHeaders:
                self.TableData[ColNameOrNumb] = DataArray
            else:
                print("Error! Column name \"%s\" was not found in the data table!"%(ColNameOrNumb))
                return
        else:
            if ColNameOrNumb >= 0 and ColNameOrNumb < len(self.TableHeaders):
                self.TableData[int(ColNameOrNumb)] = DataArray
            else:
                print("Error! Column number %d is out of allowed range (0 - %d)!"%(int(ColNameOrNumb),len(self.TableHeaders)-1))
                return
    # 
    def saveAs(self, OutputFilePath, OverWrite = False):
        if os.path.isfile(OutputFilePath):
            if OverWrite == True:
                os.system("mv %s %s"%(OutputFilePath, OutputFilePath+'.backup'))
                self.FitsStruct.writeto(OutputFilePath)
                print("Output to %s! (A backup has been created as %s)"%(OutputFilePath, OutputFilePath+'.backup'))
            else:
                print("We will not overwrite unless you specify saveAs(OverWrite=True)!")
        else:
            self.FitsStruct.writeto(OutputFilePath)
            print("Output to %s!"%(OutputFilePath))



def recognize_Filter(input_str, second_str=''):
    Filter_Name = ""
    Filter_Wave = numpy.nan # um
    if type(input_str) is str:
        # 
        # Check if AB_25, mJy things
        if input_str.endswith('_AB_25'):
            input_str = input_str.replace('_AB_25', '')
        if input_str.endswith('_mJy'):
            input_str = input_str.replace('_mJy', '')
        # 
        # Match filters
        if input_str == 'FLUX_GALEX_NUV' or input_str == 'FLUXERR_GALEX_NUV': 
            Filter_Name = 'GALEX NUV'
            Filter_Wave = 2313.9 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str.find('_U_KPNO')>=0: 
            Filter_Name = 'KPNO'
            Filter_Wave = 0.3593  # um, Skelton 2014ApJS..214...24S Table 6
        # 
        elif input_str.find('_G_Keck')>=0: 
            Filter_Name = 'Keck LRIS'
            Filter_Wave = 0.4751  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_R_Keck')>=0: 
            Filter_Name = 'Keck LRIS'
            Filter_Wave = 0.6819  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_Rs_Keck')>=0: 
            Filter_Name = 'Keck LRIS'
            Filter_Wave = 0.6819  # um, Skelton 2014ApJS..214...24S Table 6
        # 
        #Pattern_Subaru = re.compile("^([IN][AB])([0-9][0-9][0-9])[^0-9]*")
        #Matched_Subaru = Pattern_Subaru.match(input_str)
        #if Matched_Subaru: 
        #    Filter_Name = 'Subaru Suprime-Cam'
        #    Filter_Wave = float(Matched_Subaru.groups()[1]) * 1e-3 # um
        # 
        elif input_str == 'IA427_MAG_AUTO' or input_str == 'IA427_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4263.4e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA464_MAG_AUTO' or input_str == 'IA464_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4635.1e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA484_MAG_AUTO' or input_str == 'IA484_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4849.2e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA505_MAG_AUTO' or input_str == 'IA505_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5062.5e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA527_MAG_AUTO' or input_str == 'IA527_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5261.1e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA574_MAG_AUTO' or input_str == 'IA574_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5764.8e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA624_MAG_AUTO' or input_str == 'IA624_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 6233.1e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA679_MAG_AUTO' or input_str == 'IA679_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 6781.1e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA709_MAG_AUTO' or input_str == 'IA709_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7073.6e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA738_MAG_AUTO' or input_str == 'IA738_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7361.6e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA767_MAG_AUTO' or input_str == 'IA767_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7684.9e-4  # um, Laigle 2016 Table 1
        elif input_str == 'IA827_MAG_AUTO' or input_str == 'IA827_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 8244.5e-4  # um, Laigle 2016 Table 1
        elif input_str == 'NB711_MAG_AUTO' or input_str == 'NB711_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7119.9e-4  # um, Laigle 2016 Table 1
        elif input_str == 'NB816_MAG_AUTO' or input_str == 'NB816_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 8149.4e-4  # um, Laigle 2016 Table 1
        # 
        elif input_str == 'B_MAG_AUTO' or input_str == 'B_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4458.3 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'V_MAG_AUTO' or input_str == 'V_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5477.8 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'r_MAG_AUTO' or input_str == 'r_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 6288.7 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'ip_MAG_AUTO' or input_str == 'ip_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7683.9 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'zpp_MAG_AUTO' or input_str == 'zpp_MAGERR_AUTO': 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 9105.7 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str.find('_B_Subaru')>=0: 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 0.4448  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_V_Subaru')>=0: 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 0.5470  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_R_Subaru')>=0: 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 0.6276  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_i_Subaru')>=0: 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 0.7671  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_z_Subaru')>=0: 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 0.9028  # um, Skelton 2014ApJS..214...24S Table 6
        # 
        elif input_str.find('_J_Subaru')>=0: 
            Filter_Name = 'Subaru MOIRCS'
            Filter_Wave = 1.2517  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_H_Subaru')>=0: 
            Filter_Name = 'Subaru MOIRCS'
            Filter_Wave = 1.6347  # um, Skelton 2014ApJS..214...24S Table 6
        elif input_str.find('_Ks_Subaru')>=0: 
            Filter_Name = 'Subaru MOIRCS'
            Filter_Wave = 2.1577  # um, Skelton 2014ApJS..214...24S Table 6
        # 
        elif input_str == 'Y_MAG_AUTO' or input_str == 'Y_MAGERR_AUTO': 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 10214.2 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'J_MAG_AUTO' or input_str == 'J_MAGERR_AUTO': 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 12534.6 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'H_MAG_AUTO' or input_str == 'H_MAGERR_AUTO': 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 16453.4 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'Ks_MAG_AUTO' or input_str == 'Ks_MAGERR_AUTO': 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 21539.9 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str == 'Hw_MAG_AUTO' or input_str == 'Hw_MAGERR_AUTO': 
            Filter_Name = 'CFHT WIRCAM'
            Filter_Wave = 16311.4 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str == 'Ksw_MAG_AUTO' or input_str == 'Ksw_MAGERR_AUTO': 
            Filter_Name = 'CFHT WIRCAM'
            Filter_Wave = 21590.4 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str == 'u_MAG_AUTO' or input_str == 'u_MAGERR_AUTO': 
            Filter_Name = 'CFHT MegaCam'
            Filter_Wave = 3823.3 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str == 'FLUX_435W' or input_str == 'FLUXERR_435W' or input_str == 'f_f435w' or input_str == 'df_f435w' or input_str == 'e_f435w': 
            Filter_Name = 'HST ACS'
            Filter_Wave = 0.435 # um
        elif input_str == 'FLUX_606W' or input_str == 'FLUXERR_606W' or input_str == 'f_f606w' or input_str == 'df_f606w' or input_str == 'e_f606w': 
            Filter_Name = 'HST ACS'
            Filter_Wave = 0.606 # um
        elif input_str == 'FLUX_775W' or input_str == 'FLUXERR_775W' or input_str == 'f_f775w' or input_str == 'df_f775w' or input_str == 'e_f775w': 
            Filter_Name = 'HST ACS'
            Filter_Wave = 0.775 # um
        elif input_str == 'FLUX_814W' or input_str == 'FLUXERR_814W' or input_str == 'f_f814w' or input_str == 'df_f814w' or input_str == 'e_f814w': 
            Filter_Name = 'HST ACS'
            Filter_Wave = 0.814 # um
        elif input_str == 'FLUX_850W' or input_str == 'FLUXERR_850W' or input_str == 'f_f850lp' or input_str == 'df_f850lp' or input_str == 'e_f850lp': 
            Filter_Name = 'HST ACS'
            Filter_Wave = 0.850 # um
        elif input_str == 'FLUX_125W' or input_str == 'FLUXERR_125W' or input_str == 'f_f125w' or input_str == 'df_f125w' or input_str == 'e_f125w': 
            Filter_Name = 'HST WFC3'
            Filter_Wave = 1.25 # um
        elif input_str == 'FLUX_140W' or input_str == 'FLUXERR_140W' or input_str == 'f_f140w' or input_str == 'df_f140w' or input_str == 'e_f140w': 
            Filter_Name = 'HST WFC3'
            Filter_Wave = 1.40 # um
        elif input_str == 'FLUX_160W' or input_str == 'FLUXERR_160W' or input_str == 'f_f160w' or input_str == 'df_f160w' or input_str == 'e_f160w': 
            Filter_Name = 'HST WFC3'
            Filter_Wave = 1.60 # um
        # 
        elif input_str == 'f_irac1' or input_str == 'df_irac1' or input_str == 'e_irac1': 
            Filter_Name = 'Spitzer IRAC ch1'
            Filter_Wave = 35634.3 * 1e-4 # um
        elif input_str == 'f_irac2' or input_str == 'df_irac2' or input_str == 'e_irac2': 
            Filter_Name = 'Spitzer IRAC ch2'
            Filter_Wave = 45110.1 * 1e-4 # um
        elif input_str == 'f_irac3' or input_str == 'df_irac3' or input_str == 'e_irac3': 
            Filter_Name = 'Spitzer IRAC ch3'
            Filter_Wave = 57593.4 * 1e-4 # um
        elif input_str == 'f_irac4' or input_str == 'df_irac4' or input_str == 'e_irac4': 
            Filter_Name = 'Spitzer IRAC ch4'
            Filter_Wave = 79594.9 * 1e-4 # um
        # 
        elif input_str == 'SPLASH_1_FLUX' or input_str == 'SPLASH_1_FLUX_ERR': 
            Filter_Name = 'Spitzer IRAC ch1'
            Filter_Wave = 35634.3 * 1e-4 # um
        elif input_str == 'SPLASH_2_FLUX' or input_str == 'SPLASH_2_FLUX_ERR': 
            Filter_Name = 'Spitzer IRAC ch2'
            Filter_Wave = 45110.1 * 1e-4 # um
        elif input_str == 'SPLASH_3_FLUX' or input_str == 'SPLASH_3_FLUX_ERR': 
            Filter_Name = 'Spitzer IRAC ch3'
            Filter_Wave = 57593.4 * 1e-4 # um
        elif input_str == 'SPLASH_4_FLUX' or input_str == 'SPLASH_4_FLUX_ERR': 
            Filter_Name = 'Spitzer IRAC ch4'
            Filter_Wave = 79594.9 * 1e-4 # um
        # 
        elif input_str == 'FLUX_16' or input_str == 'FLUXERR_16' or input_str == 'f16' or input_str == 'df16' or input_str == 'e16': 
            Filter_Name = 'Spitzer IRS PUI'
            Filter_Wave = 16.0 # um
        # 
        elif input_str == 'FLUX_24' or input_str == 'FLUXERR_24' or input_str == 'f24' or input_str == 'df24' or input_str == 'e24': 
            Filter_Name = 'Spitzer MIPS'
            Filter_Wave = 24.0 # um
        # 
        elif input_str == 'FLUX_100' or input_str == 'FLUXERR_100' or input_str == 'f100' or input_str == 'df100' or input_str == 'e100': 
            Filter_Name = 'Herschel PACS 100'
            Filter_Wave = 100.0 # um
        elif input_str == 'FLUX_160' or input_str == 'FLUXERR_160' or input_str == 'f160' or input_str == 'df160' or input_str == 'e160': 
            Filter_Name = 'Herschel PACS 160'
            Filter_Wave = 160.0 # um
        elif input_str == 'FLUX_250' or input_str == 'FLUXERR_250' or input_str == 'f250' or input_str == 'df250' or input_str == 'e250': 
            Filter_Name = 'Herschel SPIRE 250'
            Filter_Wave = 250.0 # um
        elif input_str == 'FLUX_350' or input_str == 'FLUXERR_350' or input_str == 'f350' or input_str == 'df350' or input_str == 'e350': 
            Filter_Name = 'Herschel SPIRE 350'
            Filter_Wave = 350.0 # um
        elif input_str == 'FLUX_500' or input_str == 'FLUXERR_500' or input_str == 'f500' or input_str == 'df500' or input_str == 'e500': 
            Filter_Name = 'Herschel SPIRE 500'
            Filter_Wave = 500.0 # um
        # 
        elif input_str == 'FLUX_850' or input_str == 'FLUXERR_850' or input_str == 'f850' or input_str == 'df850' or input_str == 'e850': 
            Filter_Name = 'SCUBA-2'
            Filter_Wave = 850.0 # um
        # 
        elif input_str == 'FLUX_1160' or input_str == 'FLUXERR_1160' or input_str == 'f1160' or input_str == 'df1160' or input_str == 'e1160': 
            Filter_Name = 'AzTEC'
            Filter_Wave = 1160.0 # um
        # 
        elif input_str == 'FLUX_20cm' or input_str == 'FLUXERR_20cm' or input_str == 'f20cm' or input_str == 'df20cm' or input_str == 'e20cm': 
            Filter_Name = 'VLA'
            Filter_Wave = 2e5 # um
        # 
    # 
    return Filter_Name, Filter_Wave





##########################################
#               MAIN PROGRAM             #
##########################################

if len(sys.argv) > 1:
    
    TableFile = sys.argv[1]
    TableData = CrabFitsTable(TableFile)
    
    SOU = [] # SOURCE Table Headers
    SED = [] # MAG/FLUX/ERROR Table Headers
    
    for TableHeader in TableData.getColumnNames():
        # 
        # SOURCE
        Pattern_SOU = re.compile("[^A-Z]*(SOURCE)[^A-Z]*")
        if Pattern_SOU.match(TableHeader):
            SOU.append(TableHeader)
        Pattern_SOU = re.compile("[^A-Z]*(OBJECT)[^A-Z]*")
        if Pattern_SOU.match(TableHeader):
            SOU.append(TableHeader)
        Pattern_SOU = re.compile("[^A-Z]*(NAME)[^A-Z]*")
        if Pattern_SOU.match(TableHeader):
            SOU.append(TableHeader)
        Pattern_SOU = re.compile("[^A-Z]*(ID)[^A-Z]*")
        if Pattern_SOU.match(TableHeader):
            SOU.append(TableHeader)
        Pattern_SOU = re.compile("[^A-Z]*(SUBID)[^A-Z]*")
        if Pattern_SOU.match(TableHeader):
            SOU.append(TableHeader)
        Pattern_SOU = re.compile("[^A-Z]*(id_1)[^A-Z]*")
        if Pattern_SOU.match(TableHeader):
            SOU.append(TableHeader)
        # 
        # SED
        # 
        # MAG
        if TableHeader.upper().find('_MAG_AUTO') > 0 : 
            SED.append({'Type':'MAG AB', 'Head':TableHeader, 'Band':TableHeader.replace('_MAG_AUTO',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'AB'})
            continue
        # MAGERR
        if TableHeader.upper().find("_MAGERR_AUTO") > 0 : 
            SED.append({'Type':'MAG ERROR', 'Head':TableHeader, 'Band':TableHeader.replace('_MAGERR_AUTO',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'AB'})
            continue
        # 
        # FLUX (like FLUX_XXX or XXX_FLUX)
        if TableHeader.upper().startswith('FLUX_') or TableHeader.upper().endswith('_FLUX'): 
            if TableHeader.upper().startswith('FLUX_RADIUS'):
                pass
            #elif TableHeader.upper().startswith('FLUX_814W'):
            #    #<TODO># Laigle 2016 FLUX_814W should be magnitude AB, and the values are consistent with NB816_MAG_AUTO at similar wavelength. 
            #    SED.append({'Type':'MAG AB', 'Head':TableHeader, 'Band':TableHeader.replace('FLUX_',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'AB'})
            elif TableHeader.upper().startswith('FLUX_24') or TableHeader.upper().startswith('FLUX_GALEX') or TableHeader.upper().startswith('SPLASH_'): 
                SED.append({'Type':'FLUX uJy', 'Head':TableHeader, 'Band':TableHeader.replace('FLUX_','').replace('_FLUX',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'uJy'})
                continue
            else:
                SED.append({'Type':'FLUX mJy', 'Head':TableHeader, 'Band':TableHeader.replace('FLUX_','').replace('_FLUX',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
                continue
        # FLUXERR (like FLUXERR_XXX or XXX_FLUXERR)
        if TableHeader.upper().startswith('FLUXERR_') or TableHeader.upper().endswith('_FLUX_ERR'): 
            if TableHeader.upper().startswith('FLUXERR_RADIUS'):
                pass
            #elif TableHeader.upper().startswith('FLUXERR_814W'):
            #    #<TODO># Laigle 2016 FLUX_814W should be magnitude AB, and the values are consistent with NB816_MAG_AUTO at similar wavelength. 
            #    SED.append({'Type':'MAG ERROR', 'Head':TableHeader, 'Band':TableHeader.replace('FLUXERR_',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'AB'})
            elif TableHeader.upper().startswith('FLUXERR_24') or TableHeader.upper().startswith('FLUXERR_GALEX') or TableHeader.upper().startswith('SPLASH_'): 
                SED.append({'Type':'FLUX ERROR uJy', 'Head':TableHeader, 'Band':TableHeader.replace('FLUXERR_','').replace('_FLUXERR','').replace('_FLUX_ERR',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'uJy'})
                continue
            else:
                SED.append({'Type':'FLUX ERROR mJy', 'Head':TableHeader, 'Band':TableHeader.replace('FLUXERR_','').replace('_FLUXERR','').replace('_FLUX_ERR',''), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
                continue
        # 
        # FLUX (f24um)
        Pattern = re.compile("^(f)([0-9]+)([cum]m)*$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX mJy', 'Head':TableHeader, 'Band':Matched.group(2), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
            if Matched.group(3) is not None: SED[len(SED)-1]['Band'] = Matched.group(2)+Matched.group(3)
            continue
        # FLUXERR (df24um)
        Pattern = re.compile("^(df)([0-9]+)([cum]m)*$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX ERROR mJy', 'Head':TableHeader, 'Band':Matched.group(2), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
            if Matched.group(3) is not None: SED[len(SED)-1]['Band'] = Matched.group(2)+Matched.group(3)
            continue
        # FLUXERR (e24um)
        Pattern = re.compile("^(e)([0-9]+)([cum]m)*$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX ERROR mJy', 'Head':TableHeader, 'Band':Matched.group(2), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
            if Matched.group(3) is not None: SED[len(SED)-1]['Band'] = Matched.group(2)+Matched.group(3)
            continue
        # 
        # FLUX (f_irac1_xxx_AB_25)
        Pattern = re.compile("^(f)_([a-zA-Z0-9_]+)*(_AB_25)$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX AB 25', 'Head':TableHeader, 'Band':Matched.group(2)+Matched.group(3), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'AB 25'})
            continue
        # FLUXERR (df_irac1_xxx_AB_25)
        Pattern = re.compile("^(df)_([a-zA-Z0-9_]+)*(_AB_25)$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX ERROR AB 25', 'Head':TableHeader, 'Band':Matched.group(2)+Matched.group(3), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'AB 25'})
            continue
        # FLUXERR (e_irac1_xxx_AB_25)
        Pattern = re.compile("^(e)_([a-zA-Z0-9_]+)*(_AB_25)$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX ERROR AB 25', 'Head':TableHeader, 'Band':Matched.group(2)+Matched.group(3), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'AB 25'})
            continue
        # 
        # FLUX (f_irac1_xxx)
        Pattern = re.compile("^(f)_([a-zA-Z0-9_]+)*$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX mJy', 'Head':TableHeader, 'Band':Matched.group(2), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
            continue
        # FLUXERR (df_irac1_xxx)
        Pattern = re.compile("^(df)_([a-zA-Z0-9_]+)*$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX ERROR mJy', 'Head':TableHeader, 'Band':Matched.group(2), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
            continue
        # FLUXERR (e_irac1_xxx)
        Pattern = re.compile("^(e)_([a-zA-Z0-9_]+)*$")
        Matched = Pattern.match(TableHeader)
        if Matched:
            SED.append({'Type':'FLUX ERROR mJy', 'Head':TableHeader, 'Band':Matched.group(2), 'Wave':numpy.nan, 'Flux':numpy.nan, 'Flux Error':numpy.nan, 'WaveUnit':'um', 'FluxUnit':'mJy'})
            continue
    # 
    #for x in SED:
    #    print(x['Band'], x['Type'])
    # 
    #for x in SOU:
    #    print(x)
    
    
    if len(sys.argv) > 2:
        
        for i in range(len(sys.argv)-2):
            SourceID = sys.argv[i+2]
            print("# Getting Source by Input %s"%(SourceID))
            # 
            # check the input SourceID type (dict or str)
            if SourceID.startswith('{') and SourceID.endswith('}'):
                #print(SourceID)
                SourceID_Dict = json.loads(SourceID.replace("'",'"')) # json only accepts double quotes
                #print(SourceID_Dict)
            else:
                SourceID_Dict = {}
                for j in range(len(SOU)):
                    SourceID_Dict[SOU[j]] = SourceID
            # 
            # loop SOURCE Table Headers and match SourceID
            SourceID_Match = []
            for j in range(len(SOU)):
                if SOU[j] in SourceID_Dict:
                    #print(SOU[j])
                    SourceID_Where = numpy.argwhere(TableData.getColumn(SOU[j]).astype(str) == str(SourceID_Dict[SOU[j]])).T # do a transpose to numpy.argwhere(), see http://stackoverflow.com/questions/33747908/output-of-numpy-wherecondition-is-not-an-array-but-a-tuple-of-arrays-why
                    if len(SourceID_Where) > 0:
                        #print(SOU[j])
                        if len(SourceID_Match) > 0:
                            SourceID_Match = numpy.intersect1d(SourceID_Match, SourceID_Where)
                            #print(SourceID_Match)
                        else:
                            SourceID_Match = SourceID_Where
                            #print(SourceID_Match)
            # 
            # found Source matched by SourceID
            #print(len(SourceID_Match))
            for j in range(len(SourceID_Match)):
                #print("# -------------------------------------------------")
                print("# Found Source at Index %s"%(SourceID_Match[j]))
                # loop SOURCE Table Headers and print
                for k in range(len(SOU)):
                    if SOU[k] in SourceID_Dict:
                        print("# Getting Column %-41s Data %s"%(str(SOU[k]), TableData.getColumn(str(SOU[k]))[SourceID_Match[j]]))
                # loop SED Data Array and convert flux
                for k in range(len(SED)):
                    FilterHead = str(SED[k]['Head'])
                    FilterName, FilterWave = recognize_Filter(FilterHead)
                    FilterFlux = TableData.getColumn(FilterHead)[SourceID_Match[j]]
                    FilterType = SED[k]['Type']
                    FilterFluxUnit = SED[k]['FluxUnit']
                    FilterWaveUnit = 'um'
                    # 
                    if FilterType == 'MAG AB':
                        FilterFlux = math.pow(10, (float(FilterFlux)/(-2.5))) * 3630.7805 * 1e3 # AB magnitude to mJy, https://en.wikipedia.org/wiki/AB_magnitude
                        FilterFluxUnit = 'mJy'
                        FilterType = 'FLUX mJy'
                    # 
                    elif FilterType == 'FLUX uJy':
                        FilterFlux = FilterFlux * 1e-3 # mJy
                        FilterFluxUnit = 'mJy'
                        FilterType = 'FLUX mJy'
                    # 
                    elif FilterType == 'FLUX ERROR uJy':
                        FilterFlux = FilterFlux * 1e-3 # mJy
                        FilterFluxUnit = 'mJy'
                        FilterType = 'FLUX ERROR mJy'
                    # 
                    elif FilterType == 'FLUX AB 25':
                        FilterFlux = FilterFlux / 10**3.44 # mJy
                        FilterFluxUnit = 'mJy'
                        FilterType = 'FLUX mJy'
                        # magAB = 25.0-2.5*log10(flux) = 8.90-2.5*log10(f_Jy) = 8.90+2.5*3-2.5*log10(f_mJy) = 16.4-2.5*log10(f_mJy)
                        # so 8.60-2.5*log10(flux) = -2.5*log10(f_mJy)
                        # so 3.44-log10(flux) = -log10(f_mJy)
                        # so f_mJy = flux/10**3.44 = flux / 2754.228703
                    # 
                    elif FilterType == 'FLUX ERROR AB 25':
                        FilterFlux = FilterFlux / 10**3.44 # mJy
                        FilterFluxUnit = 'mJy'
                        FilterType = 'FLUX ERROR mJy'
                        # see notes above.
                    # 
                    SED[k]['Wave'] = FilterWave
                    SED[k]['WaveUnit'] = FilterWaveUnit
                    SED[k]['Flux'] = FilterFlux
                    SED[k]['FluxUnit'] = FilterFluxUnit
                    SED[k]['Type'] = FilterType
                    print("# Getting Column %-20s Wave %-15.6e Flux %-15.6e FluxUnit %s"%(FilterHead, FilterWave, FilterFlux, FilterFluxUnit))
                # sort SED
                SED_unsorted = SED
                SED_sorted = sorted(SED, key=lambda k: float('-inf') if math.isnan(k['Wave']) else k['Wave']) # sort function breaks in the presence of nan, see http://stackoverflow.com/questions/4240050/python-sort-function-breaks-in-the-presence-of-nan
                #for k in SED_sorted:
                #    print(k['Wave'])
                SED = SED_sorted
                # loop SED (sorted)
                for k in range(len(SED)):
                    # check if two consequent SED items are FLUX followed by FLUX ERROR, print if that is the case.
                    FilterHead = str(SED[k]['Head'])
                    FilterWave = float(SED[k]['Wave'])
                    FilterFlux = -99
                    FilterFErr = -99
                    if k >= 1:
                        #<DEBUG># print(k, SED[k-1]['Band'], SED[k]['Band'])
                        if SED[k-1]['Band'] == SED[k]['Band']:
                            # 
                            FilterBand = str(SED[k]['Band'])
                            # 
                            if SED[k]['Type'] == 'MAG ERROR':
                                FilterFlux = float(SED[k-1]['Flux']) # mJy
                                FilterFErr = float(SED[k]['Flux']) * float(SED[k-1]['Flux']) # mJy
                                FilterFluxUnit = SED[k-1]['FluxUnit']
                            elif SED[k]['Type'] == 'FLUX ERROR uJy':
                                FilterFlux = float(SED[k-1]['Flux']) # mJy
                                FilterFErr = float(SED[k]['Flux']) * 1e-3 # mJy
                                FilterFluxUnit = SED[k-1]['FluxUnit']
                            elif SED[k]['Type'] == 'FLUX ERROR mJy':
                                FilterFlux = float(SED[k-1]['Flux']) # mJy
                                FilterFErr = float(SED[k]['Flux']) # mJy
                                FilterFluxUnit = SED[k-1]['FluxUnit']
                            # print
                            #<DEBUG># print(k, SED[k-1]['Band'], SED[k]['Band'], FilterFlux, FilterFErr)
                            if FilterWave > 0 and FilterFlux > 0 and FilterFErr > 0:
                                print("                 %-20s Wave %-15.6e Flux %-15.6e FluxError %-15.6e FluxUnit %s"%(FilterBand, FilterWave, FilterFlux, FilterFErr, FilterFluxUnit))
            # 
            if len(SourceID_Match) == 0:
                print("Could not find source accroding to SourceID %s!"%(SourceID))
        
        












