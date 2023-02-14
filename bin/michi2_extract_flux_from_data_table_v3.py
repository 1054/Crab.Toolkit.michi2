#!/usr/bin/env python3
# 
# 20190125 
# 20220826 v3
# 

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(sys.argv[0])) + os.path.sep + 'lib' + os.path.sep + 'python' + os.path.sep + 'crabtable')

from CrabTable import CrabTable

import glob
import math
import numpy
import astropy
from astropy import units
from astropy.io import fits
import re
import json
import copy





#########################################
#               Functions               #
#########################################

def get_template_filter_dict():
    template_filter_dict = {'telescope_name': '', 'filter_name': '', 'wavelength_um': numpy.nan}
    return template_filter_dict

def get_all_filter_dict(catalog_name = ''):
    Filter_Dict = {}
    Filter_Dict['f435w']    = {'telescope_name': 'HST',      'filter_name': 'ACS F435W',           'wavelength_um': 0.43179 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f606w']    = {'telescope_name': 'HST',      'filter_name': 'ACS F606W',           'wavelength_um': 0.59194 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f775w']    = {'telescope_name': 'HST',      'filter_name': 'ACS F775W',           'wavelength_um': 0.76933 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f814w']    = {'telescope_name': 'HST',      'filter_name': 'ACS F814W',           'wavelength_um': 0.8057 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f850lp']   = {'telescope_name': 'HST',      'filter_name': 'ACS F850LP',          'wavelength_um': 0.90364 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f125w']    = {'telescope_name': 'HST',      'filter_name': 'WFC3 F125W',          'wavelength_um': 1.24710 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f140w']    = {'telescope_name': 'HST',      'filter_name': 'WFC3 F140W',          'wavelength_um': 1.39240 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f160w']    = {'telescope_name': 'HST',      'filter_name': 'WFC3 F160W',          'wavelength_um': 1.53960 } # Skelton 2014ApJS..214...24S Table 6
    Filter_Dict['f115w']    = {'telescope_name': 'JWST',     'filter_name': 'NIRCam F115W',        'wavelength_um': 1.154 } # https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters
    Filter_Dict['f150w']    = {'telescope_name': 'JWST',     'filter_name': 'NIRCam F150W',        'wavelength_um': 1.501 } # 
    Filter_Dict['f277w']    = {'telescope_name': 'JWST',     'filter_name': 'NIRCam F277W',        'wavelength_um': 2.786 } # 
    Filter_Dict['f444w']    = {'telescope_name': 'JWST',     'filter_name': 'NIRCam F444W',        'wavelength_um': 4.421 } # 
    Filter_Dict['f770w']    = {'telescope_name': 'JWST',     'filter_name': 'MIRI F770W',          'wavelength_um': 7.7 } # https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-imaging
    Filter_Dict['HST-F814W'] = Filter_Dict['f814w']
    Filter_Dict['CFHT-u']   = {'telescope_name': 'CFHT',     'filter_name': 'MegaCam u',           'wavelength_um':  3823.3e-4 } # Laigle+2016 Table 1
    Filter_Dict['u']        = Filter_Dict['CFHT-u']
    Filter_Dict['B']        = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam B',       'wavelength_um':  4458.3e-4 } # Laigle+2016 Table 1
    Filter_Dict['V']        = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam V',       'wavelength_um':  5477.8e-4 } # Laigle+2016 Table 1
    Filter_Dict['r']        = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam r',       'wavelength_um':  6288.7e-4 } # Laigle+2016 Table 1
    Filter_Dict['ip']       = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam i',       'wavelength_um':  7683.9e-4 } # Laigle+2016 Table 1
    Filter_Dict['zp']       = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam z',       'wavelength_um':  9105.7e-4 } # Laigle+2016 Table 1
    Filter_Dict['zpp']      = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam z',       'wavelength_um':  9105.7e-4 } # Laigle+2016 Table 1
    Filter_Dict['HSC-g']    = {'telescope_name': 'Subaru',   'filter_name': 'HSC g',               'wavelength_um':  4847e-4 } # Weaver+2021 Table 1
    Filter_Dict['HSC-r']    = {'telescope_name': 'Subaru',   'filter_name': 'HSC r',               'wavelength_um':  6219e-4 } # Weaver+2021 Table 1
    Filter_Dict['HSC-i']    = {'telescope_name': 'Subaru',   'filter_name': 'HSC i',               'wavelength_um':  7699e-4 } # Weaver+2021 Table 1
    Filter_Dict['HSC-z']    = {'telescope_name': 'Subaru',   'filter_name': 'HSC z',               'wavelength_um':  8894e-4 } # Weaver+2021 Table 1
    Filter_Dict['HSC-y']    = {'telescope_name': 'Subaru',   'filter_name': 'HSC y',               'wavelength_um':  9761e-4 } # Weaver+2021 Table 1
    Filter_Dict['yHSC']     = Filter_Dict['HSC-y']
    Filter_Dict['UVISTA-Y']     = {'telescope_name': 'VISTA',  'filter_name': 'VIRCAM Y',          'wavelength_um': 10214.2e-4 } # Laigle+2016 Table 1
    Filter_Dict['UVISTA-J']     = {'telescope_name': 'VISTA',  'filter_name': 'VIRCAM J',          'wavelength_um': 12534.6e-4 } # Laigle+2016 Table 1
    Filter_Dict['UVISTA-H']     = {'telescope_name': 'VISTA',  'filter_name': 'VIRCAM H',          'wavelength_um': 16453.4e-4 } # Laigle+2016 Table 1
    Filter_Dict['UVISTA-Ks']    = {'telescope_name': 'VISTA',  'filter_name': 'VIRCAM K',          'wavelength_um': 21539.9e-4 } # Laigle+2016 Table 1
    Filter_Dict['UVISTA-NB118'] = {'telescope_name': 'VISTA',  'filter_name': 'VIRCAM NB118',      'wavelength_um': 11909e-4 } # Weaver+2021 Table 1
    Filter_Dict['CFHT-H']       = {'telescope_name': 'CFHT',     'filter_name': 'WIRCam H',            'wavelength_um': 16453.4e-4 } # Laigle+2016 Table 1
    Filter_Dict['CFHT-Ks']      = {'telescope_name': 'CFHT',     'filter_name': 'WIRCam K',            'wavelength_um': 21539.9e-4 } # Laigle+2016 Table 1
    Filter_Dict['Y']            = Filter_Dict['UVISTA-Y']
    Filter_Dict['J']            = Filter_Dict['UVISTA-J']
    Filter_Dict['H']            = Filter_Dict['UVISTA-H']
    Filter_Dict['Ks']           = Filter_Dict['UVISTA-Ks']
    Filter_Dict['Hw']           = Filter_Dict['CFHT-H']
    Filter_Dict['K']            = Filter_Dict['CFHT-Ks']
    Filter_Dict['Ksw']          = Filter_Dict['CFHT-Ks']
    Filter_Dict['IA427']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA427',   'wavelength_um':  4263.4e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA464']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA464',   'wavelength_um':  4635.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA484']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA484',   'wavelength_um':  4849.2e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA505']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA505',   'wavelength_um':  5062.5e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA527']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA527',   'wavelength_um':  5261.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA574']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA574',   'wavelength_um':  5764.8e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA624']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA624',   'wavelength_um':  6233.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA679']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA679',   'wavelength_um':  6781.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA709']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA709',   'wavelength_um':  7073.6e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA738']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA738',   'wavelength_um':  7361.6e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA767']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA767',   'wavelength_um':  7684.9e-4 } # Laigle+2016 Table 1
    Filter_Dict['IA827']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IA827',   'wavelength_um':  8244.5e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB427']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB427',   'wavelength_um':  4263.4e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB464']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB464',   'wavelength_um':  4635.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB484']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB484',   'wavelength_um':  4849.2e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB505']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB505',   'wavelength_um':  5062.5e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB527']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB527',   'wavelength_um':  5261.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB574']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB574',   'wavelength_um':  5764.8e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB624']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB624',   'wavelength_um':  6233.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB679']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB679',   'wavelength_um':  6781.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB709']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB709',   'wavelength_um':  7073.6e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB738']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB738',   'wavelength_um':  7361.6e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB767']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB767',   'wavelength_um':  7684.9e-4 } # Laigle+2016 Table 1
    Filter_Dict['IB827']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam IB827',   'wavelength_um':  8244.5e-4 } # Laigle+2016 Table 1
    Filter_Dict['NB711']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam NB711',   'wavelength_um':  7119.9e-4 } # Laigle+2016 Table 1
    Filter_Dict['NB816']    = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam NB816',   'wavelength_um':  8149.4e-4 } # Laigle+2016 Table 1
    # 
    if catalog_name.find('3DHST'):
        Filter_Dict['u']  = {'telescope_name': 'KPNO',     'filter_name': 'u',              'wavelength_um': 0.35929 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['B']  = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam B',  'wavelength_um': 0.44480 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['G']  = {'telescope_name': 'Keck',     'filter_name': 'LRIS G',         'wavelength_um': 0.47508 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['V']  = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam V',  'wavelength_um': 0.54702 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['R']  = {'telescope_name': 'Keck',     'filter_name': 'LRIS R',         'wavelength_um': 0.62755 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['Rs'] = {'telescope_name': 'Keck',     'filter_name': 'LRIS Rs',        'wavelength_um': 0.68186 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['i']  = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam i',  'wavelength_um': 0.76712 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['z']  = {'telescope_name': 'Subaru',   'filter_name': 'Suprime Cam z',  'wavelength_um': 0.90282 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['J']  = {'telescope_name': 'Subaru',   'filter_name': 'MOIRCS J',       'wavelength_um': 1.25170 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['H']  = {'telescope_name': 'Subaru',   'filter_name': 'MOIRCS H',       'wavelength_um': 1.63470 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['Ks'] = {'telescope_name': 'Subaru',   'filter_name': 'MOIRCS Ks',      'wavelength_um': 2.15770 }  # Skelton 2014ApJS..214...24S Table 6
        Filter_Dict['b'] = Filter_Dict['B']
        Filter_Dict['g'] = Filter_Dict['G']
        Filter_Dict['v'] = Filter_Dict['V']
        Filter_Dict['r'] = Filter_Dict['R']
        Filter_Dict['rs'] = Filter_Dict['Rs']
        Filter_Dict['j'] = Filter_Dict['J']
        Filter_Dict['h'] = Filter_Dict['H']
        Filter_Dict['ks'] = Filter_Dict['Ks']
    # 
    Filter_Dict['IRAC_ch1'] = {'telescope_name': 'Spitzer',  'filter_name': 'IRAC ch1',            'wavelength_um': 35634.3e-4 } # Laigle+2016 Table 1
    Filter_Dict['IRAC_ch2'] = {'telescope_name': 'Spitzer',  'filter_name': 'IRAC ch2',            'wavelength_um': 45110.1e-4 } # Laigle+2016 Table 1
    Filter_Dict['IRAC_ch3'] = {'telescope_name': 'Spitzer',  'filter_name': 'IRAC ch3',            'wavelength_um': 57593.4e-4 } # Laigle+2016 Table 1
    Filter_Dict['IRAC_ch4'] = {'telescope_name': 'Spitzer',  'filter_name': 'IRAC ch4',            'wavelength_um': 79594.9e-4 } # Laigle+2016 Table 1
    Filter_Dict['IRAC1']    = Filter_Dict['IRAC_ch1']
    Filter_Dict['IRAC2']    = Filter_Dict['IRAC_ch2']
    Filter_Dict['IRAC3']    = Filter_Dict['IRAC_ch3']
    Filter_Dict['IRAC4']    = Filter_Dict['IRAC_ch4']
    Filter_Dict['irac1']    = Filter_Dict['IRAC_ch1']
    Filter_Dict['irac2']    = Filter_Dict['IRAC_ch2']
    Filter_Dict['irac3']    = Filter_Dict['IRAC_ch3']
    Filter_Dict['irac4']    = Filter_Dict['IRAC_ch4']
    Filter_Dict['ch1']      = Filter_Dict['IRAC_ch1']
    Filter_Dict['ch2']      = Filter_Dict['IRAC_ch2']
    Filter_Dict['ch3']      = Filter_Dict['IRAC_ch3']
    Filter_Dict['ch4']      = Filter_Dict['IRAC_ch4']
    Filter_Dict['SPLASH_1'] = Filter_Dict['IRAC_ch1']
    Filter_Dict['SPLASH_2'] = Filter_Dict['IRAC_ch2']
    Filter_Dict['SPLASH_3'] = Filter_Dict['IRAC_ch3']
    Filter_Dict['SPLASH_4'] = Filter_Dict['IRAC_ch4']
    Filter_Dict['SPLASH1']  = Filter_Dict['IRAC_ch1']
    Filter_Dict['SPLASH2']  = Filter_Dict['IRAC_ch2']
    Filter_Dict['SPLASH3']  = Filter_Dict['IRAC_ch3']
    Filter_Dict['SPLASH4']  = Filter_Dict['IRAC_ch4']
    Filter_Dict['16']       = {'telescope_name': 'Spitzer',       'filter_name': 'IRS PUI',             'wavelength_um': 16.0 }    # Liu2018, Jin2018
    Filter_Dict['24']       = {'telescope_name': 'Spitzer',       'filter_name': 'MIPS 24',             'wavelength_um': 24.0 }    # Liu2018, Jin2018
    Filter_Dict['70']       = {'telescope_name': 'Herschel',      'filter_name': 'PACS 70',             'wavelength_um': 70.0 }    # Liu2018, Jin2018
    Filter_Dict['100']      = {'telescope_name': 'Herschel',      'filter_name': 'PACS 100',            'wavelength_um': 100.0 }   # Liu2018, Jin2018
    Filter_Dict['160']      = {'telescope_name': 'Herschel',      'filter_name': 'PACS 160',            'wavelength_um': 160.0 }   # Liu2018, Jin2018
    Filter_Dict['250']      = {'telescope_name': 'Herschel',      'filter_name': 'SPIRE 250',           'wavelength_um': 250.0 }   # Liu2018, Jin2018
    Filter_Dict['350']      = {'telescope_name': 'Herschel',      'filter_name': 'SPIRE 350',           'wavelength_um': 350.0 }   # Liu2018, Jin2018
    Filter_Dict['500']      = {'telescope_name': 'Herschel',      'filter_name': 'SPIRE 500',           'wavelength_um': 500.0 }   # Liu2018, Jin2018
    Filter_Dict['850']      = {'telescope_name': 'JCMT',          'filter_name': 'SCUBA2 850',          'wavelength_um': 850.0 }   # Liu2018, Jin2018
    Filter_Dict['1100']     = {'telescope_name': 'IRAM30m',       'filter_name': 'AzTEC 1mm',           'wavelength_um': 1100.0 }  # Liu2018, Jin2018
    Filter_Dict['1160']     = {'telescope_name': '(Penner2011)',  'filter_name': 'AzTEC+MAMBO 1.16mm',  'wavelength_um': 1160.0 }  # Liu2018, Jin2018
    Filter_Dict['1200']     = {'telescope_name': 'IRAM30m',       'filter_name': 'MAMBO 1.2mm',         'wavelength_um': 1200.0 }  # Liu2018, Jin2018
    Filter_Dict['2000']     = {'telescope_name': 'IRAM30m',       'filter_name': 'GISMO 2mm',           'wavelength_um': 2000.0 }  # Liu2018, Jin2018
    Filter_Dict['3GHz']     = {'telescope_name': 'VLA',           'filter_name': '3GHz',                'wavelength_um': 1e5 }     # Liu2018, Jin2018
    Filter_Dict['1.4GHz']   = {'telescope_name': 'VLA',           'filter_name': '1.4GHz',              'wavelength_um': 2e5 }     # Liu2018, Jin2018
    Filter_Dict['10cm']     = Filter_Dict['3GHz']
    Filter_Dict['20cm']     = Filter_Dict['1.4GHz']
    # 
    # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
    Filter_Dict['sdss.up'] = {'telescope_name': 'SDSS', 'filter_name':'SDSS.u+', 'wavelength_um': 0.350102} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
    Filter_Dict['sdss.gp'] = {'telescope_name': 'SDSS', 'filter_name':'SDSS.g+', 'wavelength_um': 0.476311} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
    Filter_Dict['sdss.rp'] = {'telescope_name': 'SDSS', 'filter_name':'SDSS.r+', 'wavelength_um': 0.624698} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
    Filter_Dict['sdss.ip'] = {'telescope_name': 'SDSS', 'filter_name':'SDSS.i+', 'wavelength_um': 0.771828} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
    Filter_Dict['sdss.zp'] = {'telescope_name': 'SDSS', 'filter_name':'SDSS.z+', 'wavelength_um': 1.082983} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
    Filter_Dict['ukidss.z'] = {'telescope_name': 'UKIRT', 'filter_name':'UKIDSS.Z', 'wavelength_um': 0.8817} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=UKIRT&gname2=UKIDSS&asttype=
    Filter_Dict['ukidss.y'] = {'telescope_name': 'UKIRT', 'filter_name':'UKIDSS.Y', 'wavelength_um': 1.0305} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=UKIRT&gname2=UKIDSS&asttype=
    Filter_Dict['ukidss.j'] = {'telescope_name': 'UKIRT', 'filter_name':'UKIDSS.J', 'wavelength_um': 1.2483} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=UKIRT&gname2=UKIDSS&asttype=
    Filter_Dict['ukidss.h'] = {'telescope_name': 'UKIRT', 'filter_name':'UKIDSS.H', 'wavelength_um': 1.6313} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=UKIRT&gname2=UKIDSS&asttype=
    Filter_Dict['ukidss.k'] = {'telescope_name': 'UKIRT', 'filter_name':'UKIDSS.K', 'wavelength_um': 2.2010} # http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=UKIRT&gname2=UKIDSS&asttype=
    Filter_Dict['WISE1'] = {'telescope_name': 'WISE', 'filter_name': 'WISE1', 'wavelength_um': 3.3526} # 
    Filter_Dict['WISE2'] = {'telescope_name': 'WISE', 'filter_name': 'WISE2', 'wavelength_um': 4.6028} # 
    Filter_Dict['WISE3'] = {'telescope_name': 'WISE', 'filter_name': 'WISE3', 'wavelength_um': 11.5608} # 
    Filter_Dict['WISE4'] = {'telescope_name': 'WISE', 'filter_name': 'WISE4', 'wavelength_um': 22.0883} # 
    Filter_Dict['IRAC1'] = {'telescope_name': 'Spitzer', 'filter_name': 'IRAC1', 'wavelength_um': 3.6} # 
    Filter_Dict['IRAC2'] = {'telescope_name': 'Spitzer', 'filter_name': 'IRAC1', 'wavelength_um': 4.5} # 
    Filter_Dict['IRAC3'] = {'telescope_name': 'Spitzer', 'filter_name': 'IRAC1', 'wavelength_um': 5.8} # 
    Filter_Dict['IRAC4'] = {'telescope_name': 'Spitzer', 'filter_name': 'IRAC1', 'wavelength_um': 8.0} # 
    Filter_Dict['MIPS24']  = {'telescope_name': 'Spitzer', 'filter_name': 'MIPS24', 'wavelength_um': 24.} # 
    Filter_Dict['MIPS70']  = {'telescope_name': 'Spitzer', 'filter_name': 'MIPS70', 'wavelength_um': 70.} # 
    Filter_Dict['MIPS160'] = {'telescope_name': 'Spitzer', 'filter_name': 'MIPS160', 'wavelength_um': 160.} # 
    Filter_Dict['PACS70']  = {'telescope_name': 'Herschel', 'filter_name': 'PACS70', 'wavelength_um': 71.8} # dzliu
    Filter_Dict['PACS100'] = {'telescope_name': 'Herschel', 'filter_name': 'PACS100', 'wavelength_um': 103.0} # dzliu
    Filter_Dict['PACS160'] = {'telescope_name': 'Herschel', 'filter_name': 'PACS160', 'wavelength_um': 167.0} # dzliu
    Filter_Dict['SPIRE250'] = {'telescope_name': 'Herschel', 'filter_name': 'SPIRE250', 'wavelength_um': 252.0} # dzliu
    Filter_Dict['SPIRE350'] = {'telescope_name': 'Herschel', 'filter_name': 'SPIRE350', 'wavelength_um': 353.0} # dzliu
    Filter_Dict['SPIRE500'] = {'telescope_name': 'Herschel', 'filter_name': 'SPIRE500', 'wavelength_um': 511.0} # dzliu
    # 
    return Filter_Dict


def find_matched_filter(matched_band, matched_suffix='', special_file_name=''):
    all_filter_dict = get_all_filter_dict()
    matched_filter_dict = None
    # if the extracted band name matches one of the filter_dict
    if matched_band not in all_filter_dict and matched_band.lower() in all_filter_dict:
        matched_band = matched_band.lower()
    if matched_band not in all_filter_dict and matched_band.upper() in all_filter_dict:
        matched_band = matched_band.upper()
    if matched_band in all_filter_dict:
        matched_filter_dict = copy.copy(all_filter_dict[matched_band])
        if matched_suffix != '': 
            matched_filter_dict['filter_name'] += matched_suffix
    # elif the extracted band name is a frequency value with GHz unit
    elif re.match(r'^([0-9.]+)GHz$', matched_band):
        temp_frequency_value = float(re.sub(r'^([0-9.]+)GHz$', r'\1', matched_band))
        matched_filter_dict = get_template_filter_dict()
        matched_filter_dict['wavelength_um'] = 2.99792458e5 / temp_frequency_value
        matched_filter_dict['filter_name'] = re.sub(r'^_+',r'',matched_suffix) # '(%s)'%(re.sub(r'^_+',r'',matched_suffix))
    # elif the extracted band name is a wavelength value with (mm,um) unit
    elif re.match(r'^([0-9.]+)(mm|um|micron)$', matched_band):
        temp_wavelength_value = float(re.sub(r'^([0-9.]+)(mm|um|micron)$', r'\1', matched_band))
        temp_wavelength_unit = re.sub(r'^([0-9.]+)(mm|um|micron)$', r'\2', matched_band)
        matched_filter_dict = get_template_filter_dict()
        matched_filter_dict['wavelength_um'] = temp_wavelength_value * dict(mm=1e3, um=1., micron=1.)[temp_wavelength_unit]
        matched_filter_dict['filter_name'] = re.sub(r'^_+',r'',matched_suffix) # '(%s)'%(re.sub(r'^_+',r'',matched_suffix))
    return matched_filter_dict


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


def recognize_Col_FLUX(input_list, special_file_name='', skip_list=None):
    # output:
    # 1. recognized col name list and 
    # 2. filter_list
    recognized_list = []
    filter_list = []
    if type(input_list) is not list:
        input_list = [input_list]
    for i, input_str in enumerate(input_list):
        if skip_list is not None:
            if input_str in skip_list:
                continue
        if type(input_str) is not str:
            input_str = str(input_str)
        # 
        # we allow following format
        #   (ANYTEXT_)FLUX(_)Ks(_ANYTEXT)
        #   (ANYTEXT_)f(_)Ks(_ANYTEXT)
        #   (ANYTEXT_)S(_)Ks(_ANYTEXT)
        #  
        matched_prefix, matched_type, matched_band, matched_suffix = ['', '', '', '']
        matcher = re.match(r'^(FLUX|f|S)_([0-9a-zA-Z.-]+?)(_.*|)$', input_str, re.IGNORECASE)
        if matcher is not None:
            matched_type, matched_band, matched_suffix = matcher.groups()
        else:
            matcher = re.match(r'^(.*_|)(FLUX|f|S)_([0-9a-zA-Z.-]+?)(_.*|)$', input_str, re.IGNORECASE)
            if matcher is not None:
                matched_prefix, matched_type, matched_band, matched_suffix = matcher.groups()
            else:
                matcher = re.match(r'^(.*_|)(FLUX|f|S)_?([0-9.]+[0-9a-zA-Z.]*?|K|Ks|ch1|ch2|ch3|ch4|irac1|irac2|irac3|irac4)(_.*|)$', input_str, re.IGNORECASE)
                if matcher is not None:
                    matched_prefix, matched_type, matched_band, matched_suffix = matcher.groups()
        # 
        if matched_prefix != '' and matched_type != '' and matched_band == '' and matched_suffix == '':
            # matched_prefix usually indicates catalog but some times also indicates band, 
            # e.g., SPLASH1_FLUX, HERSCHEL_PACS_100_FLUX
            #<TODO># 
            matched_band = matched_prefix
            matched_prefix = ''
        # 
        matched_ok = True
        if matched_band != '':
            # find filter by band name
            matched_filter_dict = find_matched_filter(matched_band, matched_suffix, special_file_name)
            # then store into lists
            if matched_filter_dict is not None:
                recognized_list.append(input_str)
                filter_list.append(copy.copy(matched_filter_dict))
                print('# Recognized column "%s" as filter %s'%(input_str, matched_filter_dict))
            else:
                # if no match, then check if these are known situations
                if special_file_name.find('3DHST') >= 0 or input_str.find('3DHST') >= 0:
                    if input_str.startswith('faper_') or input_str.startswith('eaper_') or input_str.startswith('e_') or input_str.startswith('w_') or input_str.startswith('nexp_') or input_str.startswith('flux_radius') or input_str.find('_flag') >= 0:
                        continue
                # if not known situations, then report error
                print('# Error! Failed to recognize column {!r} as prefix {!r}, type {!r}, band {!r} and suffix {!r}'.format(input_str, matched_prefix, matched_type, matched_band, matched_suffix))
                matched_ok = False
        # check
        if matched_ok == False:
            raise ValueError('Error! Failed to recognize column! Please check the printed message!')
    # 
    return recognized_list, filter_list


def recognize_Col_FLUXERR(input_list, special_file_name=''):
    # output:
    # 1. recognized col name list and 
    # 2. filter_list
    recognized_list = []
    filter_list = []
    filter_dict = get_all_filter_dict()
    if type(input_list) is not list:
        input_list = [input_list]
    for input_str in input_list:
        if type(input_str) is not str:
            input_str = str(input_str)
        # 
        # we allow following format
        #   (ANYTEXT_)FLUX(_)Ks(_ANYTEXT)
        #   (ANYTEXT_)f(_)Ks(_ANYTEXT)
        #   (ANYTEXT_)S(_)Ks(_ANYTEXT)
        #  
        matched_prefix, matched_type, matched_band, matched_suffix = ['', '', '', '']
        matcher = re.match(r'^(FLUXERR|FLUX_ERR|E_FLUX|df|eS|e_S|e)_([0-9a-zA-Z.-]+?)(_.*|)$', input_str, re.IGNORECASE)
        if matcher is not None:
            matched_type, matched_band, matched_suffix = matcher.groups()
        else:
            matcher = re.match(r'^(.*_|)(FLUXERR|FLUX_ERR|E_FLUX|df|eS|e_S|e)_([0-9a-zA-Z.-]+?)(_.*|)$', input_str, re.IGNORECASE)
            if matcher is not None:
                matched_prefix, matched_type, matched_band, matched_suffix = matcher.groups()
            else:
                matcher = re.match(r'^(.*_|)(FLUXERR|FLUX_ERR|E_FLUX|df|eS|e_S|e)_?([0-9.]+[0-9a-zA-Z.]*?|K|Ks|ch1|ch2|ch3|ch4|irac1|irac2|irac3|irac4)(_.*|)$', input_str, re.IGNORECASE)
                if matcher is not None:
                    matched_prefix, matched_type, matched_band, matched_suffix = matcher.groups()
        # 
        if matched_prefix != '' and matched_type != '' and matched_band == '' and matched_suffix == '':
            # matched_prefix usually indicates catalog but some times also indicates band, 
            # e.g., SPLASH1_FLUX, HERSCHEL_PACS_100_FLUX
            #<TODO># 
            matched_band = matched_prefix
            matched_prefix = ''
        # 
        matched_ok = True
        if matched_band != '':
            # find filter by band name
            matched_filter_dict = find_matched_filter(matched_band, matched_suffix, special_file_name)
            # then store into lists
            if matched_filter_dict is not None:
                recognized_list.append(input_str)
                filter_list.append(copy.copy(matched_filter_dict))
                print('# Recognized column "%s" as filter %s'%(input_str, matched_filter_dict))
            else:
                # if no match, then check if these are known situations
                if special_file_name.find('3DHST') >= 0 or input_str.find('3DHST') >= 0:
                    if input_str.startswith('faper_') or input_str.startswith('eaper_') or input_str.startswith('e_') or input_str.startswith('w_') or input_str.startswith('nexp_') or input_str.startswith('flux_radius') or input_str.find('_flag') >= 0:
                        continue
                # if not known situations, then report error
                print('# Error! Failed to recognize column {!r} as prefix {!r}, type {!r}, band {!r} and suffix {!r}'.format(input_str, matched_prefix, matched_type, matched_band, matched_suffix))
                matched_ok = False
        # check
        if matched_ok == False:
            raise ValueError('Error! Failed to recognize column! Please check the printed message!')
    # 
    return recognized_list, filter_list












##########################################
#               MAIN PROGRAM             #
##########################################

if len(sys.argv) <= 1:
    
    print('Usage: michi2_extract_flux_from_data_table.py catalog_input.fits 23434')
    sys.exit()


# Read input args
DataFile = ''
SourceID_Inputs = []
MaxSNR = None
CatID = None # a string as a suffix '_from_cat_%s'%(CatID)
FluxUnit_Input = None
iarg = 1
while iarg <= (len(sys.argv)-1):
    reg_match = re.match(r'^[-]+(.*)$', sys.argv[iarg])
    is_option = (reg_match is not None)
    arg_str = sys.argv[iarg]
    if is_option:
        arg_str = '--'+reg_match.group(1).upper()
        if arg_str == '--MAXSNR':
            iarg += 1
            if iarg <= (len(sys.argv)-1):
                MaxSNR = float(sys.argv[iarg])
                print('# Setting max SNR limit to %s'%(MaxSNR))
        elif arg_str == '--CATID':
            iarg += 1
            if iarg <= (len(sys.argv)-1):
                CatID = sys.argv[iarg]
                print('# Setting catalog ID to %s'%(CatID)) # if set, then we will append "_from_cat_%d" to the output file name.
        elif arg_str == '--ID':
            iarg += 1
            if iarg <= (len(sys.argv)-1):
                SourceID_Inputs.append(sys.argv[iarg])
                print('# Selecting source ID %s'%(sys.argv[iarg])) # if set, then we will only select source ID matching this string.
        elif arg_str == '--FLUXUNIT':
            iarg += 1
            if iarg <= (len(sys.argv)-1):
                FluxUnit_Input = sys.argv[iarg]
                print('# Setting FluxUnit to %s'%(FluxUnit_Input)) # if set, then we will set the FluxUnit
        
    else:
        if DataFile == '':
            DataFile = sys.argv[iarg]
        else:
            SourceID_Inputs.append(sys.argv[iarg])
    # 
    iarg += 1


# Read data file
if DataFile != '':
    print('# Reading "%s"'%(DataFile))
    DataTable = CrabTable(DataFile, verbose=0)
    
    Col_Source = recognize_Col_Source(DataTable.getColumnNames(), special_file_name=DataFile)
    Col_ID = recognize_Col_ID(DataTable.getColumnNames(), special_file_name=DataFile)
    Col_FLUXERR, Flt_FLUXERR = recognize_Col_FLUXERR(DataTable.getColumnNames(), special_file_name=DataFile)
    Col_FLUX, Flt_FLUX= recognize_Col_FLUX(DataTable.getColumnNames(), special_file_name=DataFile, skip_list=Col_FLUXERR)
    
    # special treatment
    
    if 1 == 1:
        print('# Col_Source = %s'%Col_Source)
        print('# Col_ID = %s'%Col_ID)
        print('# Col_FLUX = %s'%Col_FLUX)
        print('# Col_FLUXERR = %s'%Col_FLUXERR)
    
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
    
    SourceID_Matches = [] # an array of dict
    
    if len(SourceID_Inputs) > 0:
        
        for SourceID_Input in SourceID_Inputs:
            print('# Finding source by the input name or id "%s"'%(SourceID_Input))
            # 
            # check the SourceID_Input type (dict or str)
            if SourceID_Input.startswith('{') and SourceID_Input.endswith('}'):
                #print(SourceID_Input) # something like {'ID_Laigle':12345}
                SourceID_Dict = json.loads(SourceID_Input.replace("'",'"')) # json only accepts double quotes
                #print(SourceID_Dict)
            else:
                SourceID_Dict = {}
                for icol in Col_Source:
                    SourceID_Dict[icol] = SourceID_Input
                #print(SourceID_Dict)
            # 
            # loop SOURCE Table Headers and match SourceID_Input
            SourceID_Match = []
            for icol in Col_Source:
                if icol in SourceID_Dict:
                    Search_List = DataTable.getColumn(icol).astype(str).tolist()
                    Search_Item = str(SourceID_Dict[icol])
                    SourceID_Where = numpy.argwhere(numpy.array(Search_List) == Search_Item).flatten()
                    if len(SourceID_Where) > 0:
                        if len(SourceID_Match) > 0:
                            SourceID_Match = numpy.intersect1d(SourceID_Match, SourceID_Where) # solve multiplicity <TODO>
                        else:
                            SourceID_Match = SourceID_Where
            # 
            # check whether we found any source according to the input SourceID_Input
            if len(SourceID_Match) >= 0:
                SourceID_Dict['row'] = SourceID_Match[0]
                SourceID_Dict['source id input'] = SourceID_Input
                for icol in Col_Source:
                    SourceID_Dict[icol] = DataTable.getColumn(icol)[SourceID_Dict['row']]
                SourceID_Matches.append(SourceID_Dict)
        # 
        # print error if not found
        if len(SourceID_Matches) == 0:
            print('***********')
            print('# Error! Could not find source according to the input name or id "%s"!'%(SourceID_Input))
            print('***********')
            raise ValueError('Error occurred! Please check the above printed message!')
        
    else:
        
        # if no SourceID_Input is given by the user, we will output flux for each object
        for irow in range(DataTable.getRowNumber()):
            SourceID_Dict = {}
            SourceID_Dict['row'] = irow
            SourceID_Dict['source id input'] = 'at_row_%d'%(irow+1)
            for icol in Col_Source:
                SourceID_Dict[icol] = DataTable.getColumn(icol)[SourceID_Dict['row']]
            SourceID_Matches.append(SourceID_Dict)
        
        print('# Getting all %d sources in the data table'%(len(SourceID_Matches)))
    
    print('')
    
    # 
    # found Source matched by SourceID_Input
    # 
    for j in range(len(SourceID_Matches)):
        #print("# -------------------------------------------------")
        print("# Matched source \"%s\"."%(SourceID_Matches[j]))
        
        # prepare SED data structure
        SED = []
        
        # loop SOURCE Table Headers and print
        #for k in range(len(Col_Source)):
        #    if Col_Source[k] in SourceID_Dict:
        #        print('# Matched with column "%s" row %d "%s"'%(str(Col_Source[k]), SourceID_Matches[j]['row'], DataTable.getColumn(str(Col_Source[k]))[SourceID_Matches[j]['row']]))
        
        # loop SED Data Array and convert flux
        for k in range(len(Col_FLUX)):
            FilterHead = Col_FLUX[k]
            #print('Debug', 'recognize_Filter', Col_FLUX[k], 'special_file_name='+DataFile)
            #FilterName, FilterWave = recognize_Filter(Col_FLUX[k], special_file_name=DataFile)
            FilterName = Flt_FLUX[k]['filter_name']
            FilterWave = Flt_FLUX[k]['wavelength_um']
            FilterFlux = DataTable.getColumn(Col_FLUX[k])[SourceID_Matches[j]['row']]
            FilterFluxErr = DataTable.getColumn(Col_FLUXERR[k])[SourceID_Matches[j]['row']]
            FilterType = 'FLUX'
            FilterFluxUnit = 'mJy'
            FilterWaveUnit = 'um'
            if FluxUnit_Input is not None:
                FilterFluxUnit = FluxUnit_Input
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
        out_name = 'extracted_flux'
        out_name_for_obj = 'extracted_flux_for_obj_%s'%(re.sub(r'\W+', ' ', SourceID_Matches[j]['source id input']).strip().replace(' ','_'))
        if CatID is not None:
            out_name_for_obj = out_name_for_obj + '_from_cat_%s'%(CatID)
        fout = out_name+'.txt'
        fout_obj = out_name_for_obj+'.txt'
        fout_info = out_name+'.info'
        fout_info_obj = out_name_for_obj+'.info'
        with open(fout_info,'w') as fp:
            fp.write('cat = "%s"  # the input catalog name\n'%(DataFile))
            fp.write('row = %d  # index starting from 0\n'%(SourceID_Matches[j]['row']))
            for icol in Col_Source:
                fp.write('%s = %s\n'%(icol, SourceID_Matches[j][icol]))
            fp.close()
        # 
        with open(fout,'w') as fp:
            for k in range(len(SED)):
                FilterName = str(SED[k]['Filter'])
                FilterWave = float(SED[k]['Wave'])
                FilterFlux = float(SED[k]['Flux'])
                FilterFErr = float(SED[k]['FluxErr'])
                FilterFluxUnit = str(SED[k]['FluxUnit'])
                # MaxSNR (20180111)
                if MaxSNR is not None:
                    if (FilterFErr>0) and (FilterFlux>MaxSNR*FilterFErr):
                        FilterFErr = FilterFlux/MaxSNR
                # 3DHST AB MAG ZERO POINT IS 25, WE NEED TO DO FLUX CONVERSION #<TODO><20180128># 
                if FilterName.endswith('_3DHST'):
                    FilterFlux = FilterFlux / 2754.228703 #<TODO><20180128># 
                    FilterFErr = FilterFErr / 2754.228703 #<TODO><20180128># 
                # print
                if k == 0:
                    print("# %-20s %-18s %-18s %-12s %-s"%('Wave', 'Flux', 'FluxErr', 'FluxUnit', 'FilterName'))
                    fp.write("# %-20s %-18s %-18s %-12s %-s\n"%('Wave', 'Flux', 'FluxErr', 'FluxUnit', 'FilterName'))
                #<20180110># if FilterWave > 0 and FilterFlux > 0 and FilterFErr > 0
                if numpy.isfinite(FilterFlux) and FilterFErr > 0:
                    if FilterFlux < 0.0:
                        FilterFlux = 0.0 # 20220829, if negative flux then set to zero.
                    #print("                 %-20s Wave %-15.6e Flux %-15.6e FluxError %-15.6e FluxUnit %s"%(FilterName, FilterWave, FilterFlux, FilterFErr, FilterFluxUnit))
                    print("  %-20.8g %-18.8g %-18.8g %-12s %-s"%(FilterWave, FilterFlux, FilterFErr, FilterFluxUnit, FilterName.replace(' ','_')))
                    fp.write("  %-20.8g %-18.8g %-18.8g %-12s %-s\n"%(FilterWave, FilterFlux, FilterFErr, FilterFluxUnit, FilterName.replace(' ','_')))
            fp.close()
        # 
        os.system('cp "%s" "%s"'%(fout, fout_obj))
        os.system('cp "%s" "%s"'%(fout_info, fout_info_obj))
        print('Output to "%s"!'%(fout))
        print('Output to "%s"!'%(fout_info))
        print('Output to "%s"!'%(fout_obj))
        print('Output to "%s"!'%(fout_info_obj))
        print('')
        # 
        #break












