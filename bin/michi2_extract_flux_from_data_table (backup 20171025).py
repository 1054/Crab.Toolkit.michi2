#!/usr/bin/env python2.7
# 

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')

from CrabTable import *

import glob
import math
import numpy
import astropy
from astropy import units
from astropy.io import fits
import re
import json





#########################################
#               Functions               #
#########################################

def recognize_Col_Source(input_list, special_file_name=''):
    recognized_list = []
    if type(input_list) is not list:
        input_list = [input_list]
    for input_str in input_list:
        if type(input_str) is not str:
            input_str = str(input_str)
        Pattern = re.compile("^[_]*(SOURCE)[^a-zA-Z]*", re.IGNORECASE)
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(OBJECT)[^a-zA-Z]*", re.IGNORECASE)
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(NAME)[^a-zA-Z]*", re.IGNORECASE)
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(SOURCE_[a-zA-Z]*)[^a-zA-Z]*", re.IGNORECASE)
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(OBJECT_[a-zA-Z]*)[^a-zA-Z]*", re.IGNORECASE)
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(NAME_[a-zA-Z]*)[^a-zA-Z]*", re.IGNORECASE)
        if Pattern.match(input_str):
            recognized_list.append(input_str)
    return list(set(recognized_list))


def recognize_Col_ID(input_list, special_file_name=''):
    recognized_list = []
    if type(input_list) is not list:
        input_list = [input_list]
    for input_str in input_list:
        if type(input_str) is not str:
            input_str = str(input_str)
        Pattern = re.compile("^[_]*(ID)[^a-zA-Z]*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(id)[^a-zA-Z]*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(ID_[a-zA-Z]*)[^a-zA-Z]*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(id_[a-zA-Z]*)[^a-zA-Z]*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(NUMBER)[^a-zA-Z]*", re.IGNORECASE)
        if Pattern.match(input_str):
            recognized_list.append(input_str)
    return list(set(recognized_list))


def recognize_Col_FLUX(input_list, special_file_name=''):
    recognized_list = []
    if type(input_list) is not list:
        input_list = [input_list]
    for input_str in input_list:
        if type(input_str) is not str:
            input_str = str(input_str)
        Pattern = re.compile("^[_]*(FLUX_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(f)([0-9]+).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(f_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[a-zA-Z]+(_FLUX_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(f)(ch[0-9]+)")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
    # 
    if special_file_name.find('Laigle')>=0:
        special_remove_list = ['FLUX_RADIUS', 'FLUX_814W', 'FLUX_XMM_0.5_2', 'FLUX_XMM_2_10', 'FLUX_XMM_5_10', 'FLUX_CHANDRA_0.5_2', 'FLUX_CHANDRA_2_10', 'FLUX_CHANDRA_0.5_10', 'FLUX_NUSTAR_3_24', 'FLUX_NUSTAR_3_8', 'FLUX_NUSTAR_8_24']
        # 'FLUX_814W' is actually MAG instead of FLUX in Laigle+ catalog
        for special_remove_item in special_remove_list:
            if special_remove_item in recognized_list:
                recognized_list.remove(special_remove_item)
    # 
    return recognized_list


def recognize_Col_FLUXERR(input_list, special_file_name=''):
    recognized_list = []
    if type(input_list) is not list:
        input_list = [input_list]
    for input_str in input_list:
        if type(input_str) is not str:
            input_str = str(input_str)
        Pattern = re.compile("^[_]*(FLUXERR_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(FLUX_ERR_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(df)([0-9]+).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(df_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[a-zA-Z]+(_FLUXERR_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(df)(ch[0-9]+)")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
    # 
    if special_file_name.find('Laigle')>=0:
        special_remove_list = ['FLUXERR_814W', 'FLUXERR_NUSTAR_3_24', 'FLUXERR_NUSTAR_3_8', 'FLUXERR_NUSTAR_8_24']
        # 'FLUX_814W' is actually MAG instead of FLUX in Laigle+ catalog
        for special_remove_item in special_remove_list:
            if special_remove_item in recognized_list:
                recognized_list.remove(special_remove_item)
    # 
    return recognized_list


def recognize_Col_MAG(input_list, special_file_name=''):
    recognized_list = []
    if type(input_list) is not list:
        input_list = [input_list]
    for input_str in input_list:
        if type(input_str) is not str:
            input_str = str(input_str)
        Pattern = re.compile("^[_]*(MAG_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(m)([0-9]+).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(m_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
    return recognized_list


def recognize_Col_MAGERR(input_list, special_file_name=''):
    recognized_list = []
    if type(input_list) is not list:
        input_list = [input_list]
    for input_str in input_list:
        if type(input_str) is not str:
            input_str = str(input_str)
        Pattern = re.compile("^[_]*(MAGERR_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(MAG_ERR_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(e)([0-9]+).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
        Pattern = re.compile("^[_]*(e_).*")
        if Pattern.match(input_str):
            recognized_list.append(input_str)
    return recognized_list


def recognize_Filter(input_str, special_file_name=''):
    Filter_Name = ''
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
        elif input_str.startswith('IA427_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4263.4e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA464_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4635.1e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA484_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4849.2e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA505_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5062.5e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA527_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5261.1e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA574_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5764.8e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA624_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 6233.1e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA679_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 6781.1e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA709_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7073.6e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA738_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7361.6e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA767_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7684.9e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('IA827_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 8244.5e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('NB711_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7119.9e-4  # um, Laigle 2016 Table 1
        elif input_str.startswith('NB816_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 8149.4e-4  # um, Laigle 2016 Table 1
        # 
        elif input_str.startswith('B_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 4458.3 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('V_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 5477.8 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('r_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 6288.7 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('ip_'): 
            Filter_Name = 'Subaru Suprime-Cam'
            Filter_Wave = 7683.9 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('zpp_'): 
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
        elif input_str.startswith('Y_'): 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 10214.2 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('J_'): 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 12534.6 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('H_'): 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 16453.4 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('Ks_'): 
            Filter_Name = 'VISTA VIRCAM'
            Filter_Wave = 21539.9 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str.startswith('Hw_'): 
            Filter_Name = 'CFHT WIRCAM'
            Filter_Wave = 16311.4 * 1e-4 # um, Laigle 2016 Table 1
        elif input_str.startswith('Ksw_'): 
            Filter_Name = 'CFHT WIRCAM'
            Filter_Wave = 21590.4 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str.startswith('u_'): 
            Filter_Name = 'CFHT MegaCam'
            Filter_Wave = 3823.3 * 1e-4 # um, Laigle 2016 Table 1
        # 
        elif input_str == 'f_irac1' or input_str == 'df_irac1' or input_str == 'e_irac1' or input_str == 'fch1' or input_str == 'dfch1': 
            Filter_Name = 'Spitzer IRAC ch1'
            Filter_Wave = 35634.3 * 1e-4 # um
        elif input_str == 'f_irac2' or input_str == 'df_irac2' or input_str == 'e_irac2' or input_str == 'fch2' or input_str == 'dfch2': 
            Filter_Name = 'Spitzer IRAC ch2'
            Filter_Wave = 45110.1 * 1e-4 # um
        elif input_str == 'f_irac3' or input_str == 'df_irac3' or input_str == 'e_irac3' or input_str == 'fch3' or input_str == 'dfch3': 
            Filter_Name = 'Spitzer IRAC ch3'
            Filter_Wave = 57593.4 * 1e-4 # um
        elif input_str == 'f_irac4' or input_str == 'df_irac4' or input_str == 'e_irac4' or input_str == 'fch4' or input_str == 'dfch4': 
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
        elif input_str == 'FLUX_IRAC1' or input_str == 'FLUXERR_IRAC1': 
            Filter_Name = 'Spitzer IRAC ch1'
            Filter_Wave = 35634.3 * 1e-4 # um
        elif input_str == 'FLUX_IRAC2' or input_str == 'FLUXERR_IRAC2': 
            Filter_Name = 'Spitzer IRAC ch2'
            Filter_Wave = 45110.1 * 1e-4 # um
        elif input_str == 'FLUX_IRAC3' or input_str == 'FLUXERR_IRAC3': 
            Filter_Name = 'Spitzer IRAC ch3'
            Filter_Wave = 57593.4 * 1e-4 # um
        elif input_str == 'FLUX_IRAC4' or input_str == 'FLUXERR_IRAC4': 
            Filter_Name = 'Spitzer IRAC ch4'
            Filter_Wave = 79594.9 * 1e-4 # um
        # 
        elif input_str == 'FLUX_MIPS24' or input_str == 'FLUXERR_MIPS24' or input_str == 'f24' or input_str == 'df24': 
            Filter_Name = 'Spitzer MIPS 24'
            Filter_Wave = 24.0 # um
        # 
        elif input_str == 'FLUX_24' or input_str == 'FLUXERR_24' or input_str == 'f24' or input_str == 'df24': 
            Filter_Name = 'Spitzer MIPS 24'
            Filter_Wave = 24.0 # um
        # 
        elif input_str == 'FLUX_16' or input_str == 'FLUXERR_16' or input_str == 'f16' or input_str == 'df16': 
            Filter_Name = 'Spitzer IRS PUI 16'
            Filter_Wave = 16.0 # um
        # 
        elif input_str == 'FLUX_K' or input_str == 'FLUXERR_K' or input_str == 'fK' or input_str == 'dfK': 
            Filter_Name = 'unknown K band'
            Filter_Wave = 2.15 # um
        # 
        elif input_str == 'FLUX_814W' or input_str == 'FLUXERR_814W': 
            Filter_Name = 'HST ACS F814W'
            Filter_Wave = 814.0e-4 # um
        # 
        # 
        if numpy.isnan(Filter_Wave):
            Pattern_FLUX_1 = re.compile("(FLUX[ERR]*_)([0-9Ee.+-]+)([^0-9Ee.+-]*.*)")
            Pattern_FLUXERR_1 = re.compile("(FLUXERR_)([0-9Ee.+-]+)([^0-9Ee.+-]*.*)")
            Pattern_FLUX_2 = re.compile("[^A-Z]*(f)([0-9Ee.+-]+)([^0-9Ee.+-]*)")
            Pattern_FLUXERR_2 = re.compile("[^A-Z]*(df)([0-9Ee.+-]+)([^0-9Ee.+-]*)")
            Matched_FLUX_1 = Pattern_FLUX_1.match(input_str)
            Matched_FLUXERR_1 = Pattern_FLUXERR_1.match(input_str)
            Matched_FLUX_2 = Pattern_FLUX_2.match(input_str)
            Matched_FLUXERR_2 = Pattern_FLUXERR_2.match(input_str)
            Matched = None
            if Matched is None:
                if Matched_FLUX_1:
                    Matched = Matched_FLUX_1
            if Matched is None:
                if Matched_FLUXERR_1:
                    Matched = Matched_FLUXERR_1
            if Matched is None:
                if Matched_FLUX_2:
                    Matched = Matched_FLUX_2
            if Matched is None:
                if Matched_FLUXERR_2:
                    Matched = Matched_FLUXERR_2
            if Matched: 
                Filter_Name = 'unknown'
                Filter_Wave = float(Matched.group(2))
                # try to convert Filter_Wave unit
                if Matched.group(3) is not None: 
                    if Matched.group(3).strip().startswith('cm'):
                        Filter_Wave = Filter_Wave * 1e4 # convert from cm to um
                    elif Matched.group(3).strip().startswith('GHz'):
                        Filter_Wave = 2.99792458e5 / Filter_Wave # convert from GHz to um
                    elif Matched.group(3) == 'W':
                        Filter_Wave = Filter_Wave / 1e4 # convert from AA to um
                # try to guess Filter_Name
                if Matched.group(3) is not None: 
                    if Matched.group(3).upper().find('ALMA')>=0:
                        Filter_Name = 'ALMA'
                    elif Matched.group(3).upper().find('SCUBA2')>=0:
                        Filter_Name = 'JCMT SCUBA2'
                    elif Matched.group(3).upper().find('SCUBA2')>=0 and Matched.group(3).upper().find('AZTEC')>=0:
                        Filter_Name = 'AzTEC SCUBA2'
                    elif Matched.group(2) == '1.4' and Matched.group(3).startswith('GHz'):
                        Filter_Name = 'VLA 1.4GHz'
                    elif Matched.group(2) == '3' and Matched.group(3).startswith('GHz'):
                        Filter_Name = 'VLA 3GHz'
                    elif Matched.group(2) == '10' and Matched.group(3).startswith('cm'):
                        Filter_Name = 'VLA 3GHz'
                if Filter_Name == 'unknown':
                    if Matched.group(2) =='70':
                        Filter_Name = 'Herschel PACS 70'
                    elif Matched.group(2) =='100':
                        Filter_Name = 'Herschel PACS 100'
                    elif Matched.group(2) =='160':
                        Filter_Name = 'Herschel PACS 160'
                    elif Matched.group(2) =='250':
                        Filter_Name = 'Herschel SPIRE 250'
                    elif Matched.group(2) =='350':
                        Filter_Name = 'Herschel SPIRE 350'
                    elif Matched.group(2) =='500':
                        Filter_Name = 'Herschel SPIRE 500'
                    elif Matched.group(2) =='850':
                        Filter_Name = 'JCMT SCUBA2 850'
    # 
    if special_file_name.find('Laigle')>=0:
        Filter_Name = Filter_Name + ' (Laigle)'
    # 
    return Filter_Name, Filter_Wave









##########################################
#               MAIN PROGRAM             #
##########################################

if len(sys.argv) == 1:
    
    print('Usage: michi2_extract_flux_from_data_table.py catalog_input.fits 23434')
    sys.exit()

else:

    DataFile = sys.argv[1]
    print('# Reading "%s"'%(DataFile))
    DataTable = CrabTable(DataFile, verbose=0)
    
    Col_Source  = recognize_Col_Source(DataTable.getColumnNames(), special_file_name=DataFile)
    Col_ID      = recognize_Col_ID(DataTable.getColumnNames(), special_file_name=DataFile)
    Col_FLUX    = recognize_Col_FLUX(DataTable.getColumnNames(), special_file_name=DataFile)
    Col_FLUXERR = recognize_Col_FLUXERR(DataTable.getColumnNames(), special_file_name=DataFile)
    
    # special treatment
    
    if 1 == 0:
        print('Col_Source = %s'%Col_Source)
        print('Col_ID = %s'%Col_ID)
        print('Col_FLUX = %s'%Col_FLUX)
        print('Col_FLUXERR = %s'%Col_FLUXERR)
    
    if Col_Source:
        if Col_ID:
            #<BUGGY># Col_Source = Col_Source.extend(Col_ID)
            Col_Source.extend(Col_ID)
    else:
        if Col_ID:
            Col_Source = Col_ID
    
    if Col_Source is None:
        print('******************************************************************************')
        print('Error! Could not determine Source or ID columns in the input catalog!')
        print('Column names: %s'%(DataTable.getColumnNames()))
        print('Col_Source = %s'%Col_Source)
        print('Col_ID = %s'%Col_ID)
        print('Col_FLUX = %s'%Col_FLUX)
        print('Col_FLUXERR = %s'%Col_FLUXERR)
        sys.exit()
    
    if len(Col_FLUX) != len(Col_FLUXERR):
        print('******************************************************************************')
        print('Error! FLUX columns and FLUXERR columns do not match in dimension! Please check the input catalog!')
        for k in range(max(len(Col_FLUX),len(Col_FLUXERR))):
            if k < len(Col_FLUX) and k < len(Col_FLUXERR):
                print('        %-30s %-30s'%(Col_FLUX[k], Col_FLUXERR[k]))
            elif k < len(Col_FLUX):
                print('        %-30s %-30s'%(Col_FLUX[k], ' '))
            elif k < len(Col_FLUXERR):
                print('        %-30s %-30s'%(' ', Col_FLUXERR[k]))
        sys.exit()
    
    SourceID_Match = []
    
    if len(sys.argv) > 2:
        
        for i in range(len(sys.argv)-2):
            SourceID_Input = sys.argv[i+2]
            print('# Getting Source by the input name or id "%s"'%(SourceID_Input))
            # 
            # check the SourceID_Input type (dict or str)
            if SourceID_Input.startswith('{') and SourceID_Input.endswith('}'):
                #print(SourceID_Input) # something like {'ID_Laigle':12345}
                SourceID_Dict = json.loads(SourceID_Input.replace("'",'"')) # json only accepts double quotes
                #print(SourceID_Dict)
            else:
                SourceID_Dict = {}
                for j in range(len(Col_Source)):
                    SourceID_Dict[Col_Source[j]] = SourceID_Input
                #print(SourceID_Dict)
            # 
            # loop SOURCE Table Headers and match SourceID_Input
            SourceID_Match = []
            for j in range(len(Col_Source)):
                if Col_Source[j] in SourceID_Dict:
                    #print(Col_Source[j])
                    SourceID_Where = numpy.argwhere(DataTable.getColumn(Col_Source[j]).astype(str) == str(SourceID_Dict[Col_Source[j]])).T # do a transpose to numpy.argwhere(), see http://stackoverflow.com/questions/33747908/output-of-numpy-wherecondition-is-not-an-array-but-a-tuple-of-arrays-why
                    if len(SourceID_Where) > 0:
                        SourceID_Where = SourceID_Where[0]
                        if len(SourceID_Match) > 0:
                            SourceID_Match = numpy.intersect1d(SourceID_Match, SourceID_Where)
                        else:
                            SourceID_Match = SourceID_Where
                    #print SourceID_Where
                    #print SourceID_Match
        # 
        # check whether we found any source according to the input SourceID_Input
        if len(SourceID_Match) == 0:
            print('# Warning! Could not find source according to the input name or id "%s"!'%(SourceID_Input))
        
    else:
        
        # if no SourceID_Input is given by the user, we will output flux for each object
        SourceID_Match = range(DataTable.getRowNumber())
        SourceID_Dict = {}
    
    # 
    # found Source matched by SourceID_Input
    # 
    #print(len(SourceID_Match))
    for j in range(len(SourceID_Match)):
        #print("# -------------------------------------------------")
        print("# Found Source at row number %s (starting from 1)."%(SourceID_Match[j]+1))
        
        # prepare SED data structure
        SED = []
        
        # loop SOURCE Table Headers and print
        for k in range(len(Col_Source)):
            if Col_Source[k] in SourceID_Dict:
                print('# Matched with Column "%s" "%s"'%(str(Col_Source[k]), DataTable.getColumn(str(Col_Source[k]))[SourceID_Match[j]]))
        
        # loop SED Data Array and convert flux
        for k in range(len(Col_FLUX)):
            FilterHead = Col_FLUX[k]
            FilterName, FilterWave = recognize_Filter(Col_FLUX[k], special_file_name=DataFile)
            FilterFlux = DataTable.getColumn(Col_FLUX[k])[SourceID_Match[j]]
            FilterFluxErr = DataTable.getColumn(Col_FLUXERR[k])[SourceID_Match[j]]
            FilterType = 'FLUX'
            FilterFluxUnit = 'mJy'
            FilterWaveUnit = 'um'
            # 
            TempColumn = DataTable.TableColumns[DataTable.getColumnIndex(Col_FLUX[k])]
            if type(TempColumn) is astropy.io.fits.column.Column:
                if type(TempColumn.unit) is str:
                    #print(TempColumn.unit)
                    if TempColumn.unit != '':
                        FilterFluxUnit = TempColumn.unit
            # 
            if FilterFluxUnit == 'uJy':
                FilterFlux = FilterFlux / 1e3
                FilterFluxErr = FilterFluxErr / 1e3
                FilterFluxUnit = 'mJy'
            
            # 
            #if FilterType == 'MAG AB':
            #    FilterFlux = math.pow(10, (float(FilterFlux)/(-2.5))) * 3630.7805 * 1e3 # AB magnitude to mJy, https://en.wikipedia.org/wiki/AB_magnitude
            #    FilterFluxUnit = 'mJy'
            #    FilterType = 'FLUX mJy'
            #elif FilterType == 'FLUX uJy':
            #    FilterFlux = FilterFlux * 1e-3 # mJy
            #    FilterFluxUnit = 'mJy'
            #    FilterType = 'FLUX mJy'
            #
            #elif FilterType == 'FLUX AB 25':
            #    FilterFlux = FilterFlux / 10**3.44 # mJy
            #    FilterFluxUnit = 'mJy'
            #    FilterType = 'FLUX mJy'
            #    # magAB = 25.0-2.5*log10(flux) = 8.90-2.5*log10(f_Jy) = 8.90+2.5*3-2.5*log10(f_mJy) = 16.4-2.5*log10(f_mJy)
            #    # so 8.60-2.5*log10(flux) = -2.5*log10(f_mJy)
            #    # so 3.44-log10(flux) = -log10(f_mJy)
            #    # so f_mJy = flux/10**3.44 = flux / 2754.228703
            # 
            # 
            SED_k = {}
            SED_k['Head'] = FilterHead
            SED_k['Filter'] = FilterName
            SED_k['Wave'] = FilterWave
            SED_k['WaveUnit'] = FilterWaveUnit
            SED_k['Flux'] = FilterFlux
            SED_k['FluxErr'] = FilterFluxErr
            SED_k['FluxUnit'] = FilterFluxUnit
            SED_k['Type'] = FilterType
            SED.append(SED_k)
            #print("# Getting Column %-20s Wave %-15.6e Flux %-15.6e FluxUnit %s"%(FilterHead, FilterWave, FilterFlux, FilterFluxUnit))
        
        # sort SED
        SED_unsorted = SED
        SED_sorted = sorted(SED, key=lambda k: float('-inf') if math.isnan(k['Wave']) else k['Wave']) # sort function breaks in the presence of nan, see http://stackoverflow.com/questions/4240050/python-sort-function-breaks-in-the-presence-of-nan
        SED = SED_sorted
        
        # loop SED (sorted) and output to file
        fout = 'extracted_flux.txt'
        if len(SourceID_Match)>1:
            #fout = 'extracted_flux_at_row_%s.txt'%(SourceID_Match[j]+1)
            fout = 'extracted_flux_for_obj_%d.txt'%(j)
        fp = open(fout,'w')
        for k in range(len(SED)):
            FilterName = str(SED[k]['Filter'])
            FilterWave = float(SED[k]['Wave'])
            FilterFlux = float(SED[k]['Flux'])
            FilterFErr = float(SED[k]['FluxErr'])
            FilterFluxUnit = str(SED[k]['FluxUnit'])
            # print
            if k == 0:
                print("# %-20s %-18s %-18s %-12s %-s"%('Wave', 'Flux', 'FluxErr', 'FluxUnit', 'FilterName'))
                fp.write("# %-20s %-18s %-18s %-12s %-s\n"%('Wave', 'Flux', 'FluxErr', 'FluxUnit', 'FilterName'))
            if FilterWave > 0 and FilterFlux > 0 and FilterFErr > 0:
                #print("                 %-20s Wave %-15.6e Flux %-15.6e FluxError %-15.6e FluxUnit %s"%(FilterName, FilterWave, FilterFlux, FilterFErr, FilterFluxUnit))
                print("  %-20.8g %-18.8g %-18.8g %-12s %-s"%(FilterWave, FilterFlux, FilterFErr, FilterFluxUnit, FilterName.replace(' ','_')))
                fp.write("  %-20.8g %-18.8g %-18.8g %-12s %-s\n"%(FilterWave, FilterFlux, FilterFErr, FilterFluxUnit, FilterName.replace(' ','_')))
        fp.close()
        print('Output to "%s"!'%(fout))
        # 
        #break












