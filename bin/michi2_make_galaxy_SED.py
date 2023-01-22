#!/usr/bin/env python
# 
# Aim: 
#      This code template SEDs and combine them into a mock panchromatic galaxy SED
#      then outputs wavelength and flux to a two-column text file. 
#      
# Inputs:
#      See Usage
# 
# 

import os
import sys
import numpy
import astropy
import zipfile # zipfile.ZipFile
import tempfile
import shutil
from contextlib import contextmanager
import astropy.io.ascii as asciitable
from copy import copy

Script_dir = os.path.dirname(os.path.abspath(__file__))

Lib_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'lib_SED')
if not os.path.isdir(Lib_dir): 
    print('Error! Lib_dir "%s" does not exist! Please make sure you have completely downloaded this code from "github.com/1054/Crab.Toolkit.michi2"!'%(Lib_dir))
    sys.exit()

sys.path.append(Script_dir)
from michi2_read_lib_SEDs import lib_file_get_header, lib_file_get_data_block, lib_file_get_data_block_quick, check_array, nearest_interp, spline, integrate_vLv



####################################
#               FUNC               #
####################################

@contextmanager
def get_temp_dir():
    # get temporary dir
    # see https://stackoverflow.com/questions/21257782/how-to-remove-file-when-program-exits
    path = tempfile.mkdtemp()
    try:
        yield path
    finally:
        try:
            shutil.rmtree(path)
        except IOError:
            sys.stderr.write('Failed to clean up temp dir {}'.format(path))



def find_SED_lib_by_given_parameters(Lib_file_path, Lib_params_dict, verbose = True): 
    # 
    if verbose: print('find_SED_lib_by_given_parameters ...')
    if verbose: print('Lib_file_path: ', Lib_file_path)
    if verbose: print('Lib_params_dict: ', Lib_params_dict)
    # 
    # Get temp dir for storing lib files
    with get_temp_dir() as Temp_dir:
        # 
        # check the Lib_file_path, if it is a zip file then extract it to tmp directory
        if not os.path.isfile(Lib_file_path) and os.path.isfile(Lib_file_path+'.zip'):
            with zipfile.ZipFile(Lib_file_path+'.zip', 'r') as zipObj:
                Lib_file_name = os.path.basename(Lib_file_path)
                zipObj.extract(Lib_file_name, Temp_dir)
                Lib_file_path = os.path.join(Temp_dir, Lib_file_name)
                if verbose: print('Lib_file_path: ', Lib_file_path)
        # 
        # Get lib header
        Lib_header = lib_file_get_header(Lib_file_path)
        if verbose: print('Lib_header: ', Lib_header)
        # 
        # Find key column
        Lib_keys = Lib_params_dict
        Lib_key_cols = {} # the column index of each key (param)
        for Lib_key in Lib_keys:
            iLibHeaderColumn = -1
            while iLibHeaderColumn < 10:
                #<TODO># only support 10 LIB parameters
                if 'TPAR%d'%(iLibHeaderColumn+1) in Lib_header:
                    if Lib_header['TPAR%d'%(iLibHeaderColumn+1)] == Lib_key:
                        Lib_key_cols[Lib_key] = iLibHeaderColumn
                        break
                iLibHeaderColumn = iLibHeaderColumn + 1
            if iLibHeaderColumn == 10:
                print('Error! Could not find "%s" in SED lib "%s"!'%(Lib_key, Lib_file_path))
                sys.exit()
                return [], [], {}
        # 
        # Prepare to list SED lib rows for each key (param)
        Lib_rows = []
        Lib_key_rows = {} # the data block first row index of each key (param)
        Lib_key_vals = {} # the data block first row value of each key (param)
        Lib_key_all_rows = {} # all the rows indices of each key (param)
        Lib_key_all_vals = {} # all the rows values of each key (param)
        for Lib_key in Lib_keys:
            Lib_key_rows[Lib_key] = []
            Lib_key_vals[Lib_key] = []
            Lib_key_all_rows[Lib_key] = []
            Lib_key_all_vals[Lib_key] = []
        # 
        # List SED lib rows by matching the key value to the given value
        #for Lib_key in Lib_keys:
        #    Lib_key_all_vals[Lib_key] = os.pipe('CrabTableReadColumn "%s" %d'%(Lib_file_path, 2+Lib_key_cols[Lib_key]))
        #    Lib_key_all_rows[Lib_key] = numpy.arange(Lib_key_all_vals[Lib_key])
        # List the first line of each data block (template) to get the parameters of each template
        Lib_data_blocks = lib_file_get_data_block_quick(Lib_file_path, 0, Lib_header = Lib_header, every_data_line = int(Lib_header['NVAR1']), max_data_line_count = int(Lib_header['NVAR2']), verbose = False)
        #asciitable.write(Lib_data_blocks[0:10,:], sys.stdout, Writer=asciitable.FixedWidth)
        iLibDataLine = 0
        if verbose: sys.stdout.write('Reading Lib file: 0%')
        if verbose: sys.stdout.flush()
        while iLibDataLine < int(Lib_header['NVAR1']) * int(Lib_header['NVAR2']):
            #Lib_data_block = lib_file_get_data_block_quick(Lib_file_path, iLibDataLine, Lib_header = Lib_header, max_data_line_count = 1, verbose = False)
            Lib_data_block = Lib_data_blocks[int(iLibDataLine/int(Lib_header['NVAR1']))]
            Lib_rows.append(iLibDataLine)
            for Lib_key in Lib_keys:
                Lib_key_rows[Lib_key].append(iLibDataLine) # this is actually the same as Lib_rows, just duplicated for each key (param)
                Lib_key_vals[Lib_key].append(Lib_data_block[2+Lib_key_cols[Lib_key]]) # here '2+' is because there are two wavelength and flux columns
            iLibDataLine = iLibDataLine + int(Lib_header['NVAR1'])
            #if verbose:
            #    if (iLibDataLine / int(Lib_header['NVAR1'])) % int(float(Lib_header['NVAR2']) / 300.0) == 0:
            #        sys.stdout.write(' %.2f%%'%(float(iLibDataLine) / (float(Lib_header['NVAR1']) * float(Lib_header['NVAR2']) / 100.0)))
            #        sys.stdout.flush()
        if verbose: sys.stdout.write(' 100%\n')
        if verbose: sys.stdout.flush()
        # 
        # Print
        if verbose: print('Lib_keys: ', Lib_keys)
        if verbose: print('Lib_rows: ', Lib_rows)
        if verbose: print('Lib_key_cols: ', Lib_key_cols)
        if verbose: print('Lib_key_rows: ', Lib_key_rows)
        if verbose: print('Lib_key_vals: ', Lib_key_vals)
        # 
        # Match to Lib_keys
        Mask_lower = []
        Mask_upper = []
        for Lib_key in Lib_keys:
            if len(Mask_lower) == 0:
                Mask_lower = ((numpy.array(Lib_key_vals[Lib_key])-float(Lib_keys[Lib_key])) <= 0)
            else:
                Mask_lower = Mask_lower & ((numpy.array(Lib_key_vals[Lib_key])-float(Lib_keys[Lib_key])) <= 0)
            if len(Mask_upper) == 0:
                Mask_upper = ((numpy.array(Lib_key_vals[Lib_key])-float(Lib_keys[Lib_key])) >= 0)
            else:
                Mask_upper = Mask_upper & ((numpy.array(Lib_key_vals[Lib_key])-float(Lib_keys[Lib_key])) >= 0)
        # 
        if len(numpy.argwhere(Mask_lower)) == 0:
            print('Error! The given parameters are not in the allow range (out of the lower range)!')
            for Lib_key in Lib_keys:
                print('The input key %s = %s, allowed range is %s to %s.'%(Lib_key, Lib_keys[Lib_key], numpy.min(numpy.array(Lib_key_vals[Lib_key])), numpy.max(numpy.array(Lib_key_vals[Lib_key]))))
            sys.exit()
        if len(numpy.argwhere(Mask_upper)) == 0:
            print('Error! The given parameters are not in the allow range (out of the upper range)!')
            for Lib_key in Lib_keys:
                print('The input key %s = %s, allowed range is %s to %s.'%(Lib_key, Lib_keys[Lib_key], numpy.min(numpy.array(Lib_key_vals[Lib_key])), numpy.max(numpy.array(Lib_key_vals[Lib_key]))))
            sys.exit()
        # 
        # We must assume that SED lib parameters are monochromatically increasing
        if len(numpy.argwhere(Mask_lower)) > 0 and len(numpy.argwhere(Mask_upper)) > 0:
            iLibKey_lower = (numpy.argwhere(Mask_lower).flatten().tolist())[-1]
            iLibKey_upper = (numpy.argwhere(Mask_upper).flatten().tolist())[0]
            if iLibKey_lower > iLibKey_upper:
                # This happens when the Lib_key is found multiple times, and the above code will lead to 
                # iLibKey_lower takes the last matched one while
                # iLibKey_upper takes the first matched one, 
                # so we need to swap them for multiple matches case
                iLibKey_lower, iLibKey_upper = (iLibKey_upper, iLibKey_lower)
            iLibDataLine_lower = Lib_rows[iLibKey_lower]
            iLibDataLine_upper = Lib_rows[iLibKey_upper]
            if verbose: print('iLibKey_lower = ', iLibKey_lower)
            if verbose: print('iLibKey_upper = ', iLibKey_upper)
            if iLibKey_lower == iLibKey_upper:
                # Lib_data_block
                Lib_data_block = lib_file_get_data_block(Lib_file_path, iLibDataLine_lower, Lib_header = Lib_header, verbose = False)
                if verbose: print('Lib_data_block[0] = ', Lib_data_block[0])
                # No need to interpolate
                # SED_params_dict
                SED_params_dict = {}
                iLibHeaderColumn = -1
                while iLibHeaderColumn < 10:
                    #<TODO># only support 10 LIB parameters
                    if 'TPAR%d'%(iLibHeaderColumn+1) in Lib_header:
                        SED_params_dict[Lib_header['TPAR%d'%(iLibHeaderColumn+1)]] = Lib_data_block[0,2+iLibHeaderColumn] # here '2+' is because there are two wavelength and flux columns
                        #<TODO># Additional feature: if the Lib parameter column is like 'flux_*', then we store the whole array
                        if Lib_header['TPAR%d'%(iLibHeaderColumn+1)].lower().startswith('flux_'):
                            SED_params_dict[Lib_header['TPAR%d'%(iLibHeaderColumn+1)]] = Lib_data_block[:,2+iLibHeaderColumn]
                    iLibHeaderColumn = iLibHeaderColumn + 1
                if verbose: print('SED_params_dict = ', SED_params_dict)
                Wavelength_um = Lib_data_block[:,0]
                Flux_density_mJy = Lib_data_block[:,1]
                # return 
                return Wavelength_um, Flux_density_mJy, SED_params_dict
            else:
                # Lib_data_block (read two records so as to do interpolation)
                Lib_data_block_lower = lib_file_get_data_block(Lib_file_path, iLibDataLine_lower, Lib_header = Lib_header, verbose = False)
                Lib_data_block_upper = lib_file_get_data_block(Lib_file_path, iLibDataLine_upper, Lib_header = Lib_header, verbose = False)
                if verbose: print('Lib_data_block_lower[0] = ', Lib_data_block_lower[0])
                if verbose: print('Lib_data_block_upper[0] = ', Lib_data_block_upper[0])
                # do interpolation
                Multiply_factor_lower = 0.0
                Multiply_factor_upper = 0.0
                Multiply_factor_total = 0.0
                for Lib_key in Lib_keys:
                    Multiply_factor_lower = Multiply_factor_lower + numpy.abs(float(Lib_key_vals[Lib_key][iLibKey_upper]) - float(Lib_keys[Lib_key]))
                    Multiply_factor_upper = Multiply_factor_upper + numpy.abs(float(Lib_key_vals[Lib_key][iLibKey_lower]) - float(Lib_keys[Lib_key]))
                    Multiply_factor_total = Multiply_factor_total + numpy.abs(float(Lib_key_vals[Lib_key][iLibKey_upper]) - float(Lib_key_vals[Lib_key][iLibKey_lower]))
                Multiply_factor_lower = Multiply_factor_lower / Multiply_factor_total
                Multiply_factor_upper = Multiply_factor_upper / Multiply_factor_total
                if verbose: print('Multiply_factor_lower = ', Multiply_factor_lower)
                if verbose: print('Multiply_factor_upper = ', Multiply_factor_upper)
                Wavelength_um = Lib_data_block_lower[:,0]
                Flux_density_mJy = Multiply_factor_lower * Lib_data_block_lower[:,1] + Multiply_factor_upper * Lib_data_block_upper[:,1]
                # SED_params_dict
                SED_params_dict = {}
                iLibHeaderColumn = -1
                while iLibHeaderColumn < 10:
                    #<TODO># only support 10 LIB parameters
                    if 'TPAR%d'%(iLibHeaderColumn+1) in Lib_header:
                        if Lib_header['TPAR%d'%(iLibHeaderColumn+1)] in Lib_keys:
                            SED_params_dict[Lib_header['TPAR%d'%(iLibHeaderColumn+1)]] = Lib_keys[Lib_header['TPAR%d'%(iLibHeaderColumn+1)]] # if the user has given the Lib key (param) then use it. 
                        else:
                            #<TODO># Additional feature: if the Lib parameter column is like 'flux_*', then we store the whole array
                            if Lib_header['TPAR%d'%(iLibHeaderColumn+1)].lower().startswith('flux_'):
                                SED_params_dict[Lib_header['TPAR%d'%(iLibHeaderColumn+1)]] = Multiply_factor_lower * Lib_data_block_lower[:,2+iLibHeaderColumn] + Multiply_factor_upper * Lib_data_block_upper[:,2+iLibHeaderColumn] # here '2+' is because there are two wavelength and flux columns
                            else:
                                # Otherwise we store the parameter value in the first line of each SED template data block, assuming they are a single value for each SED template
                                SED_params_dict[Lib_header['TPAR%d'%(iLibHeaderColumn+1)]] = Multiply_factor_lower * Lib_data_block_lower[0,2+iLibHeaderColumn] + Multiply_factor_upper * Lib_data_block_upper[0,2+iLibHeaderColumn] # here '2+' is because there are two wavelength and flux columns
                    iLibHeaderColumn = iLibHeaderColumn + 1
                if verbose: print('SED_params_dict = ', SED_params_dict)
                # return 
                return Wavelength_um, Flux_density_mJy, SED_params_dict
    # 
    # if error, return empty arrays
    return [], [], {}











def michi2_make_galaxy_SED(\
    z:float = numpy.nan, 
    dL:float = numpy.nan, 
    Mstar:float = numpy.nan, 
    tau:float = numpy.nan, 
    Age:float = numpy.nan, 
    EBV:float = numpy.nan, 
    SFR:float = numpy.nan, 
    LIR:float = numpy.nan, 
    fAGN:float = numpy.nan, 
    TAGN:float = numpy.nan, 
    fPDR:float = numpy.nan, 
    Umin:float = numpy.nan, 
    Umean:float = numpy.nan, 
    qPAH:float = numpy.nan, 
    qIR:float = numpy.nan, 
    Wave_step:float = 0.01, # dex
    Output_file:str = '', 
    Silent = False, 
    ):
    
    # 
    # parse
    if numpy.isnan(LIR) and not numpy.isnan(SFR):
        LIR = SFR * 1e10 # Chabrier2003 IMF
    
    # 
    # print
    if Silent == False:
        print('##########################')
        print('# michi2_make_galaxy_SED #')
        print('##########################')
        print('z = %s'%(z))
        print('Mstar = %0.6e [Msun] (Chabrier IMF)'%(Mstar))
        print('Age = %0.6e [Gyr]'%(Age/1e9))
        print('E(B-V) = %0.3f [dex]'%(EBV))
        print('SFR = %s [Msun/yr] (Chabrier IMF, =LIR/1e10)'%(SFR))
        print('LIR = %0.6e [Lsun] (IR, rest-frame 8-1000um)'%(LIR))
    
    # 
    # check
    if numpy.isnan(z) or (numpy.isnan(Mstar) and numpy.isnan(LIR)):
        raise ValueError('Error! Please check the about inputs!')
        sys.exit()
    
    if numpy.isnan(EBV):
        EBV = 0.0
    if numpy.isnan(fAGN):
        fAGN = 0.0 # no AGN contribution
    if numpy.isnan(TAGN):
        TAGN = 1 # AGN type 1
    if numpy.isnan(fPDR):
        fPDR = 0.04 # 10**-1.4
    if ~numpy.isnan(Umin):
        Umin_cold = Umin
        Umin_warm = Umin
        Umean = (1-fPDR) * Umin_cold + fPDR * Umin_warm * (numpy.log(1e6/Umin_warm)/(1-Umin_warm/1e6))
    else:
        if numpy.isnan(Umean):
            #Umean = (1.0+z)**1.15 # Magdis2012 Umean versus z evolution
            Umean = 3.0 * (1.0+z)**1.80 # Bethermin2015 Umean versus z evolution
        # inversed calculation from Umean to Umin.
        Umax = 1e6
        Temp_Grid_Umin = numpy.arange(0.01,50.0+0.01,0.01)
        Temp_Grid_Umean = (1.0-fPDR) * Temp_Grid_Umin + fPDR * (numpy.log(Umax/Temp_Grid_Umin) / (1.0/Temp_Grid_Umin-1.0/Umax)) # see Draine 2007 SINGS Eq.(17)
        Umin = Temp_Grid_Umin[(numpy.argsort(numpy.abs(Temp_Grid_Umean-Umean)).flatten().tolist())[0]]
        # 
        # Magdis 2012 z~2 BzK 
        # u = (1.0+z)**1.15   for z = 0.0,3.0,0.1
        # set z = {0.00 0.04 0.04 0.30 0.30 0.65 0.65 1.00 1.00 1.30 1.30 1.75 1.75 2.25 2.25 3.00}
        # set u = {2.20 2.20 3.30 3.30 4.90 4.90 6.10 6.10 9.70 9.70 12.1 12.1 14.5 14.5 18.0 18.0}
        # 
        # Bethermin 2015 SED z~4
        # u = 3.0 * (1.0+z)**1.80    for z = 2.0,4.0,0.1
        # 
        # 
        #<TODO># currently we are limited by DL07 templates which have a max Umin of 25.0. <TODO> ask Bruce Draine for new templates (see emails with Laure Ciesla)
        #if Umin > 25.0:
        #    Umin = 25.0
    if numpy.isnan(qPAH):
        #qPAH = 0.47 # qPAH set to constant
        if z <= 1.5:
            qPAH = 0.47 # Magdis2012 BzK paper
        else:
            qPAH = 0.025 # Magdis2012 BzK paper
    if numpy.isnan(qIR):
        qIR = 2.35*(1+(z))**(-0.12) + numpy.log10(1.91) # qIR evolution from Delhaize+2017
        # 
        # De Haize et al. 2017
        # qIR = 2.35*(1+(z))**(-0.12)+lg(1.91) # changing from default qIR = 2.5 to qIR=2.35*(1+z)**(-0.12)+lg(1.91) (Magnelli 2015A%26A...573A..45M)
        # 
    
    if Silent == False:
        print('fAGN = %s'%(fAGN))
        print('TAGN = %s'%(TAGN))
        print('fPDR = %s'%(fPDR))
        print('Umean = %s'%(Umean))
        print('Umin = %s'%(Umin))
        print('qPAH = %s'%(qPAH))
        print('qIR = %s'%(qIR))
        print('Output_file = "%s"'%(Output_file))
        print('##########################')
    
    if Output_file != '':
        if not os.path.isdir(os.path.dirname(Output_file)):
            if os.path.dirname(Output_file) != '':
                if Silent == False:
                    print('os.makedirs("%s")'%(os.path.dirname(Output_file)))
                os.makedirs(os.path.dirname(Output_file))
    
    
    
    
    
    
    ####################################
    #               MAIN               #
    ####################################
    
    # Lumdist
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as Units
    import astropy.constants as Constants
    cosmo = FlatLambdaCDM(H0=69.8, Om0=0.27) # Freedman 2019 (bibcode:2019ApJ...882...34F)
    if numpy.isnan(dL):
        dL = cosmo.luminosity_distance(z).to(Units.Mpc).value # Mpc
    pi = numpy.pi
    if Silent == False:
        print('')
        print('dL = ', dL)
        print('')
    
    
    # 
    # Get stellar SED lib header
    # 
    a_stellar = 0.0
    w_stellar = None
    f_stellar = None
    if ~numpy.isnan(Mstar):
        # 
        # if no Age is given, we use BC03 constant SFH single-age snapshotlibrary
        if numpy.isnan(Age):
            Lib_stellar = os.path.join(Lib_dir, 'lib.BC03.Padova1994.BaSeL.Z0.0190.ConstSFH.Age200Myr.EBV.SED')
            w_stellar, f_stellar, p_stellar = find_SED_lib_by_given_parameters(Lib_stellar, {'EBV':EBV}, verbose = False)
        # 
        # if Age is given, we use BC03 constant SFH multi-age snapshot library
        elif numpy.isnan(tau) or numpy.isclose(tau, 0.0):
            Lib_stellar = os.path.join(Lib_dir, 'lib.BC03.Padova1994.BaSeL.Z0.0190.ConstSFH.MultiAge.EBV.SED')
            w_stellar, f_stellar, p_stellar = find_SED_lib_by_given_parameters(Lib_stellar, {'EBV':EBV, 'Age':Age/1e9}, verbose = False)
        # 
        # if both Age is given, we use FSPS tau-decline SFH library (<TODO>: we only have tau=1Gyr model...)
        else:
            if numpy.isclose(tau, 0.1):
                Lib_stellar = os.path.join(Lib_dir, 'lib.FSPS.CSP.tau.0p1Gyr.Padova.BaSeL.Z0.0190.EBV.SED')
            else:
                Lib_stellar = os.path.join(Lib_dir, 'lib.FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.EBV.SED')
            w_stellar, f_stellar, p_stellar = find_SED_lib_by_given_parameters(Lib_stellar, {'EBV':EBV, 'Age':Age/1e9}, verbose = False)
        # 
        # if Age is given, we use Philipp's BC03 tau-decline SFH library
        #else:
        #    Lib_stellar = '/Users/dzliu/Softwares/BC03/Output_SED_LIB/BC03_Constant_SFH_from_plang/lib.BC03.CstSFH.Z0.0190.EBV.SED' # 'BC03_Constant_SFH_from_plang/lib.BC03.CstSFH.Z0.0190.EBV.SED'
        #    #Lib_stellar = '/Users/dzliu/Softwares/BC03/Output_SED_LIB/BC03_Constant_SFH_various_ages/lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED' # 'BC03_Constant_SFH_various_ages/lib.BC03.2016.Padova1994.BaSeL.Z0.0190.EBV.SED'
        #    w_stellar, f_stellar, p_stellar = find_SED_lib_by_given_parameters(Lib_stellar, {'EBV':EBV, 'Age':Age/1e9}, verbose = True)
        
        if 'Mass' in p_stellar:
            # Stellar_mass = a / (3.839e33*1e26/(4*pi*dL**2*9.52140e48)) * float(p_stellar['Mass']) / (1+Redshift)
            a_stellar = Mstar * (3.839e33*1e26/(4*pi*dL**2*9.52140e48)) / float(p_stellar['Mass']) * (1+z)
            f_stellar = f_stellar * a_stellar
            w_stellar = w_stellar * (1+z)
            if Silent == False:
                print('a_stellar = %s'%(a_stellar))
                print('p_stellar[\'Mass\'] = %0.6e'%(p_stellar['Mass']))
                print('Mstar = %0.6e'%(Mstar))
        else:
            if Silent == False:
                print('Mstar = 0?')
    else:
        if Silent == False:
            print('Mstar = 0')
    
    
    # 
    # Get AGN SED lib header
    # 
    a_AGN = 0.0
    w_AGN = None
    f_AGN = None
    if fAGN > 0.0:
        Lib_AGN = os.path.join(Lib_dir, 'lib.MullaneyAGN.SED')
        w_AGN, f_AGN, p_AGN = find_SED_lib_by_given_parameters(Lib_AGN, {'AGN_TYPE':TAGN}, verbose = False)
        a_AGN =  LIR * fAGN / (p_AGN['Lbol']*4*pi*dL**2/(1+z))
        f_AGN = f_AGN * a_AGN
        w_AGN = w_AGN * (1+z)
        if Silent == False:
            print('fAGN = %0.3f%% [LIR]'%(fAGN*100.0))
            print('a_AGN = %s'%(a_AGN))
            print('L_AGN = %0.6e [Lsolar]'%(integrate_vLv(w_AGN, f_AGN, z)))
    else:
        if Silent == False:
            print('fAGN = 0% [LIR]')
    
    
    # 
    # Get warm and cold dust SED lib header
    # 
    Lib_warm_dust = os.path.join(Lib_dir, 'lib.DL07.2010.03.18.spec.2.0.HiExCom.SED')
    Lib_cold_dust = os.path.join(Lib_dir, 'lib.DL07.2010.03.18.spec.2.0.LoExCom.SED')
    w_warm_dust, f_warm_dust, p_warm_dust = find_SED_lib_by_given_parameters(Lib_warm_dust, {'Umin':Umin, 'qPAH':qPAH}, verbose = False)
    w_cold_dust, f_cold_dust, p_cold_dust = find_SED_lib_by_given_parameters(Lib_cold_dust, {'Umin':Umin, 'qPAH':qPAH}, verbose = False)
    # 
    # since
    #     Mdust_warm = a_warm_dust * dL**2 / (1+z) # Mdust #NOTE# no need to multiply a '4*pi'!
    #     Mdust_cold = a_cold_dust * dL**2 / (1+z) # Mdust #NOTE# no need to multiply a '4*pi'!
    #     fPDR = Mdust_warm / (Mdust_warm + Mdust_cold)
    # so: 
    #     1/fPDR = 1 + Mdust_cold/Mdust_warm = 1 + a_cold_dust/a_warm_dust
    # so: 
    #     a_cold_dust = a_warm_dust * (1/fPDR-1)
    # 
    # and 
    #     LIR_warm_dust = a_warm_dust * numpy.power(10,float(p_warm_dust['lgLTIR'])) * 4*pi*dL**2 / (1+z) # 40.31970 converts mJy*GHz to Lsun*Mpc-2
    #     LIR_cold_dust = a_cold_dust * numpy.power(10,float(p_cold_dust['lgLTIR'])) * 4*pi*dL**2 / (1+z) # 40.31970 converts mJy*GHz to Lsun*Mpc-2
    #     LIR_warm_dust + LIR_cold_dust = LIR
    # so:
    #     a_warm_dust * numpy.power(10,float(p_warm_dust['lgLTIR'])) + 
    #     a_cold_dust * numpy.power(10,float(p_cold_dust['lgLTIR'])) = LIR / (4*pi*dL**2 / (1+z))
    # 
    # so: 
    #     a_warm_dust * ( numpy.power(10,float(p_warm_dust['lgLTIR'])) + 
    #        (1/fPDR-1) * numpy.power(10,float(p_cold_dust['lgLTIR'])) ) = LIR / (4*pi*dL**2 / (1+z))
    # 
    # so: 
    a_warm_dust = LIR / (4*pi*dL**2 / (1+z)) / ( numpy.power(10,float(p_warm_dust['lgLTIR'])) + \
                                    (1/fPDR-1) * numpy.power(10,float(p_cold_dust['lgLTIR'])) )
    a_cold_dust = a_warm_dust * (1/fPDR-1)
    M_warm_dust = a_warm_dust * dL**2 / (1+z)
    M_cold_dust = a_cold_dust * dL**2 / (1+z)
    L_warm_dust = a_warm_dust * numpy.power(10,float(p_warm_dust['lgLTIR'])) * 4*pi*dL**2 / (1+z)
    L_cold_dust = a_cold_dust * numpy.power(10,float(p_cold_dust['lgLTIR'])) * 4*pi*dL**2 / (1+z)
    if Silent == False:
        print('a_warm_dust = %s'%(a_warm_dust))
        print('a_cold_dust = %s'%(a_cold_dust))
        print('M_warm_dust = %0.6e [Msolar]'%(M_warm_dust))
        print('M_cold_dust = %0.6e [Msolar]'%(M_cold_dust))
        print('L_warm_dust = %0.6e [Lsolar]'%(L_warm_dust))
        print('L_cold_dust = %0.6e [Lsolar]'%(L_cold_dust))
    f_warm_dust = f_warm_dust * a_warm_dust
    f_cold_dust = f_cold_dust * a_cold_dust
    w_warm_dust = w_warm_dust * (1+z)
    w_cold_dust = w_cold_dust * (1+z)
    
    
    # 
    # Get radio SED lib header
    # 
    Lib_radio = os.path.join(Lib_dir, 'lib.RadioPowerlaw.Single.SED')
    w_radio, f_radio, p_radio = find_SED_lib_by_given_parameters(Lib_radio, {'POWER_INDEX':0.8}, verbose = False) # The radio SED are normalized to f20cm(rest-frame)=1mJy
    w_radio_20cm = 20.0 # cm, rest-frame
    w_radio_1p4GHz = 29.9792458/1.4 # cm, rest-frame
    #f_radio_1p4GHz = LIR/(4*pi*dL**2)/(1+z)*1e3/3750/10**qIR # mJy
    #print('f_radio_1p4GHz = %s [mJy]'%(f_radio_1p4GHz))
    f_radio_1p4GHz = ( a_warm_dust * numpy.power(10,float(p_warm_dust['lgLTIR']))*40.31970 + \
                       a_cold_dust * numpy.power(10,float(p_cold_dust['lgLTIR']))*40.31970 )/3750/10**qIR # mJy # 40.31970 converts mJy*GHz to Lsun*Mpc-2
    if Silent == False:
        print('f_radio_1p4GHz = %s [mJy]'%(f_radio_1p4GHz))
    f_radio_20cm = f_radio_1p4GHz*(w_radio_20cm/w_radio_1p4GHz)**0.8 # mJy
    a_radio = f_radio_20cm
    f_radio = f_radio * a_radio
    w_radio = w_radio * (1+z)
    if Silent == False:
        print('f_radio_20cm = %s [mJy]'%(f_radio_20cm))
        print('a_radio = %s'%(a_radio))
    # 
    # f_radio_1p4GHz/[mJy] = LIR/[Lsolar]/3750/10**qIR
    # 
    
    
    # 
    # Combine SEDs
    # 
    log_w_SED = numpy.arange(-1.0, 5.5, Wave_step)
    w_SED = numpy.power(10,log_w_SED)
    f_SED = w_SED * 0.0
    if a_stellar > 0:
        f_SED = f_SED + spline(w_stellar, f_stellar, w_SED, xlog=1, ylog=1, fill=0.0)
    if a_AGN > 0:
        f_SED = f_SED + spline(w_AGN, f_AGN, w_SED, xlog=1, ylog=1, fill=0.0)
    if a_warm_dust > 0:
        f_SED = f_SED + spline(w_warm_dust, f_warm_dust, w_SED, xlog=1, ylog=1, fill=0.0)
    if a_cold_dust > 0:
        f_SED = f_SED + spline(w_cold_dust, f_cold_dust, w_SED, xlog=1, ylog=1, fill=0.0)
    if a_radio > 0:
        f_SED = f_SED + spline(w_radio, f_radio, w_SED, xlog=1, ylog=1, fill=0.0)
    
    
    # 
    # Output dict
    # 
    Output_dict = {}
    Output_dict['w_SED'] = w_SED
    Output_dict['f_SED'] = f_SED
    Output_dict['w_stellar'] = w_stellar
    Output_dict['f_stellar'] = f_stellar
    Output_dict['w_AGN'] = w_AGN
    Output_dict['f_AGN'] = f_AGN
    Output_dict['w_warm_dust'] = w_warm_dust
    Output_dict['f_warm_dust'] = f_warm_dust
    Output_dict['w_cold_dust'] = w_cold_dust
    Output_dict['f_cold_dust'] = f_cold_dust
    Output_dict['w_radio'] = w_radio
    Output_dict['f_radio'] = f_radio
    
    
    # 
    # Output file
    # 
    if Output_file != '':
        print('')
        asciitable.write(numpy.column_stack((w_SED,f_SED)), 
                            Output_file, 
                            Writer=asciitable.FixedWidth, 
                            names=['Wavelength_um', 'Flux_density_mJy'], 
                            overwrite=True, delimiter='  ', bookend=True)
        os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(Output_file))
        os.system('rm "%s.bak"'%(Output_file))
        print('Output to "%s"'%(Output_file))
        # 
        # Output each SED component
        # 
        Output_name, Output_extension = os.path.splitext(Output_file)
        
        if a_stellar > 0:
            if 'flux_unattenu' in p_stellar:
                asciitable.write(numpy.column_stack((w_stellar,f_stellar,p_stellar['flux_unattenu']*a_stellar)), 
                                Output_name+'_stellar'+'.txt', 
                                Writer=asciitable.FixedWidth, 
                                names=['Wavelength_um', 'Flux_density_mJy', 'Flux_density_unattenuated_mJy'], 
                                overwrite=True, delimiter='  ', bookend=True)
                print('integrate_vLv stellar SED attenuated = %0.6e [Lsolar]'%(integrate_vLv(w_stellar, f_stellar, z)))
                print('integrate_vLv stellar SED unattenuated = %0.6e [Lsolar]'%(integrate_vLv(w_stellar, p_stellar['flux_unattenu']*a_stellar, z)))
            else:
                asciitable.write(numpy.column_stack((w_stellar,f_stellar)), 
                                Output_name+'_stellar'+'.txt', 
                                Writer=asciitable.FixedWidth, 
                                names=['Wavelength_um', 'Flux_density_mJy'], 
                                overwrite=True, delimiter='  ', bookend=True)
            
            os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(Output_name+'_stellar'+'.txt'))
            os.system('rm "%s.bak"'%(Output_name+'_stellar'+'.txt'))
            print('Output to "%s"'%(Output_name+'_stellar'+'.txt'))
        
        
        if a_warm_dust > 0:
            asciitable.write(numpy.column_stack((w_warm_dust,f_warm_dust)), 
                                Output_name+'_warm_dust'+'.txt', 
                                Writer=asciitable.FixedWidth, 
                                names=['Wavelength_um', 'Flux_density_mJy'], 
                                overwrite=True, delimiter='  ', bookend=True)
            os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(Output_name+'_warm_dust'+'.txt'))
            os.system('rm "%s.bak"'%(Output_name+'_warm_dust'+'.txt'))
            print('Output to "%s"'%(Output_name+'_warm_dust'+'.txt'))
        
        
        if a_cold_dust > 0:
            asciitable.write(numpy.column_stack((w_cold_dust,f_cold_dust)), 
                                Output_name+'_cold_dust'+'.txt', 
                                Writer=asciitable.FixedWidth, 
                                names=['Wavelength_um', 'Flux_density_mJy'], 
                                overwrite=True, delimiter='  ', bookend=True)
            os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(Output_name+'_cold_dust'+'.txt'))
            os.system('rm "%s.bak"'%(Output_name+'_cold_dust'+'.txt'))
            print('Output to "%s"'%(Output_name+'_cold_dust'+'.txt'))
        
        
        if a_AGN > 0:
            asciitable.write(numpy.column_stack((w_AGN,f_AGN)), 
                                Output_name+'_AGN'+'.txt', 
                                Writer=asciitable.FixedWidth, 
                                names=['Wavelength_um', 'Flux_density_mJy'], 
                                overwrite=True, delimiter='  ', bookend=True)
            os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(Output_name+'_AGN'+'.txt'))
            os.system('rm "%s.bak"'%(Output_name+'_AGN'+'.txt'))
            print('Output to "%s"'%(Output_name+'_AGN'+'.txt'))
        
        
        if a_stellar > 0:
            asciitable.write(numpy.column_stack((w_stellar,f_stellar)), 
                                Output_name+'_stellar'+'.txt', 
                                Writer=asciitable.FixedWidth, 
                                names=['Wavelength_um', 'Flux_density_mJy'], 
                                overwrite=True, delimiter='  ', bookend=True)
            os.system('sed -i.bak -e "1s/^ /#/" "%s"'%(Output_name+'_stellar'+'.txt'))
            os.system('rm "%s.bak"'%(Output_name+'_stellar'+'.txt'))
            print('Output to "%s"'%(Output_name+'_stellar'+'.txt'))
    # 
    # Done
    # 
    return Output_dict












##################################
#              MAIN              #
##################################

if __name__ == '__main__':
    
    ####################################
    #              SYS.ARG             #
    ####################################
    
    if len(sys.argv) <= 1:
        print('Usage: ')
        print('    michi2_make_galaxy_SED.py -z        NN.N \\')
        print('                              -logMstar NN.N \\')
        print('                              -EBV      NN.N \\')
        print('                              -SFR      NN.N \\')
        print('                             #-logLIR   NN.N \\')
        print('                              -fAGN     NN.N \\')
        print('                              -TAGN     NN.N \\')
        print('                              -fPDR     NN.N \\')
        print('                             #-Umin     NN.N \\')
        print('                              -Umean    NN.N \\')
        print('                              -qPAH     NN.N \\')
        print('                              -qIR      NN.N \\')
        print('                              -Out      SSSS.txt')
        print('')
        sys.exit()
    
    z = numpy.nan
    Mstar = numpy.nan
    tau = numpy.nan
    Age = numpy.nan # yr
    EBV = numpy.nan
    SFR = numpy.nan
    LIR = numpy.nan
    fAGN = numpy.nan
    TAGN = numpy.nan
    fPDR = numpy.nan
    Umin = numpy.nan
    Umean = numpy.nan
    qPAH = numpy.nan
    qIR = numpy.nan
    Output_file = ''
    
    iarg = 1
    while iarg < len(sys.argv):
        if sys.argv[iarg].lower() == '-z':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                z = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-mstar':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                Mstar = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-logmstar':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                Mstar = numpy.power(10,float(sys.argv[iarg]))
        elif sys.argv[iarg].lower() == '-tau':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                tau = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-age':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                if sys.argv[iarg].lower().endswith('gyr'):
                    Age = float(sys.argv[iarg].lower().replace('gyr',''))*1e9
                elif sys.argv[iarg].lower().endswith('myr'):
                    Age = float(sys.argv[iarg].lower().replace('myr',''))*1e6
                elif sys.argv[iarg].lower().endswith('kyr'):
                    Age = float(sys.argv[iarg].lower().replace('kyr',''))*1e3
                elif sys.argv[iarg].lower().endswith('yr'):
                    Age = float(sys.argv[iarg].lower().replace('yr',''))
                else:
                    Age = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-logage':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                Age = numpy.power(10,float(sys.argv[iarg]))
        elif sys.argv[iarg].lower() == '-ebv':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                EBV = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-a_v':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                EBV = float(sys.argv[iarg])*0.25 # E(B-V) = 0.25 * A_V
        elif sys.argv[iarg].lower() == '-sfr':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                SFR = float(sys.argv[iarg])
                LIR = SFR * 1e10
        elif sys.argv[iarg].lower() == '-logsfr':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                SFR = numpy.power(10,float(sys.argv[iarg]))
                LIR = SFR * 1e10
        elif sys.argv[iarg].lower() == '-lir':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                LIR = float(sys.argv[iarg])
                SFR = LIR / 1e10
        elif sys.argv[iarg].lower() == '-loglir':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                LIR = numpy.power(10,float(sys.argv[iarg]))
                SFR = LIR / 1e10
        elif sys.argv[iarg].lower() == '-fagn':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                fAGN = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-tagn':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                TAGN = int(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-fpdr':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                fPDR = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-umin':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                Umin = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-umean':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                Umean = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-qpah':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                qPAH = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-qir':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                qIR = float(sys.argv[iarg])
        elif sys.argv[iarg].lower() == '-out':
            if iarg+1 < len(sys.argv):
                iarg = iarg + 1
                Output_file = str(sys.argv[iarg])
                if not Output_file.endswith('.txt'):
                    Output_file = Output_file+'.txt'
        elif sys.argv[iarg].lower() == '-test':
            z = 1.5
            Mstar = numpy.power(10,10.8)
            EBV = 0.5
            Age = 200e6 # 200 Myr
            LIR = numpy.power(10,12.8)
            SFR = LIR / 1e10
            qPAH = 0.50
            Output_file = 'test.txt'
        iarg = iarg + 1
        
        
    # 
    # 
    michi2_make_galaxy_SED(\
        z = z, 
        Mstar = Mstar, 
        tau = tau, 
        Age = Age, 
        EBV = EBV, 
        SFR = SFR, 
        LIR = LIR, 
        fAGN = fAGN, 
        TAGN = TAGN, 
        fPDR = fPDR, 
        Umin = Umin, 
        Umean = Umean, 
        qPAH = qPAH, 
        qIR = qIR, 
        Output_file = Output_file, 
    )
        
















