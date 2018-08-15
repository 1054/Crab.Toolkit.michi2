#!/usr/bin/env python3
# 
# This program will read MAGPHYS SED fitting results:
#   fit_name.sed
#   fit_name.fit
# And output result text format tables:
#   
# 
# 



# 
# Import packages
# 
import os, sys, re, json, astropy, time, itertools
import numpy as np
import astropy.io.ascii as asciitable
from astropy.table import Table
#from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)
from scipy import optimize, interpolate
from collections import Counter # Counter counts the number of occurrences of each item in a list. See -- https://stackoverflow.com/questions/30650474/python-rename-duplicates-in-list-with-progressive-numbers-without-sorting-list



# 
# Subroutine of reading magphys user_obs file
# 
def read_magphys_user_obs_file(file_path):
    return asciitable.read(file_path)



# 
# Subroutine of reading magphys filters file
# 
def read_magphys_filters_file(file_path):
    return asciitable.read(file_path)



# 
# Function for solving duplicates in a list
# 
def rename_duplicates_in_a_string_list(input_list):
    # solve duplicates in read_cols
    dupl_counts = Counter(input_list)
    for dupl_item, dupl_count in dupl_counts.items():
        if dupl_count > 1:
            for dupl_suffix in range(1, dupl_count+1): # suffix starts at 1 and increases by 1 each time
                if dupl_suffix > 1:
                    input_list[input_list.index(dupl_item)] = dupl_item + '_%d' % (dupl_suffix) # replace each appearance of s
                    # see -- https://stackoverflow.com/questions/30650474/python-rename-duplicates-in-list-with-progressive-numbers-without-sorting-list
    return input_list



# 
# Subroutine of reading magphys sed file
# 
def read_magphys_sed_file(file_path):
    fit_sed_data = {}
    line = ''
    next_line = ''
    read_cols = []
    read_vals = []
    regex_replace = re.compile(r'(\.){2,}') # replace multiple consecutive repeating character '.' # in magphys parameter header line there are lots of ... which will be replaced by white space
    with open(file_path, 'r') as fp:
        fp.seek(0,2) # go to the file end.
        EOF = fp.tell() # get the end of file location
        fp.seek(0,0) # go back to file beginning
        while fp.tell() != EOF:
            line = fp.readline()
            line = line.strip()
            if line.startswith('#'):
                # check magphys sed flux unit
                if line.find('Spectral Energy Distribution [lg(L_lambda/LoA^-1)]') >= 0:
                    fit_sed_data['__energy_unit__'] = 'lg(L_lambda/LoA^-1)'
                elif line.find('Spectral Energy Distribution [lg(Fnu/Jy)]') >= 0:
                    fit_sed_data['__energy_unit__'] = 'lg(Fnu/Jy)'
                # this line is a commented header line
                line = regex_replace.sub(' ', line[1:]).strip() # line = line[1:].replace('.',' ').strip()
                line_split = line.split()
                if len(line_split) > 0:
                    read_cols = rename_duplicates_in_a_string_list(line_split)
                else:
                    continue
            else:
                # this line is a non-commented data line
                line_split = line.split()
                if len(read_cols) > 0 and len(read_cols) == len(line_split):
                    read_vals = line_split
                    for read_col,read_val in itertools.zip_longest(read_cols,read_vals):
                        if read_col in fit_sed_data:
                            fit_sed_data[read_col].append(float(read_val))
                        else:
                            fit_sed_data[read_col] = [float(read_val)]
                        #print('fit_sed_data[\'%s\'] = %s' % (read_col, read_val)) # debug
                else:
                    continue
    if not ('__energy_unit__' in fit_sed_data):
        print('Error! Could not find "Spectral Energy Distribution [.*]" in the input MAGPHYS SED file "%d"!' % (file_path))
        sys.exit()
    # 
    return fit_sed_data



# 
# Subroutine of reading magphys fit file
# 
def read_magphys_fit_file(file_path):
    fit_fit_data = {}
    line = ''
    next_line = ''
    read_cols = []
    read_vals = []
    param_name = ''
    column_mode = True
    histogram_mode = False
    percentile_mode = False
    regex_replace = re.compile(r'(\.){2,}') # replace multiple consecutive repeating character '.' # in magphys parameter header line there are lots of ... which will be replaced by white space
    with open(file_path, 'r') as fp:
        fp.seek(0,2) # go to the file end.
        EOF = fp.tell() # get the end of file location
        fp.seek(0,0) # go back to file beginning
        while fp.tell() != EOF:
            line = fp.readline()
            line = line.strip()
            if line.startswith('#'):
                # this line is a commented header line
                line = regex_replace.sub(' ', line[1:]).strip()
                if line.find('BEST FIT MODEL: (i_sfh, i_ir, chi2, redshift)') >= 0:
                    line = 'i_sfh i_ir chi2 redshift'
                if line.find('MARGINAL PDF HISTOGRAMS FOR EACH PARAMETER') >= 0:
                    column_mode = False
                    histogram_mode = True
                if column_mode:
                    # column mode
                    line_split = line.split()
                    if len(line_split) > 0:
                        read_cols = rename_duplicates_in_a_string_list(line_split)
                    else:
                        continue
                else:
                    # histogram mode, column number is always two, but parameter name is read from current commented line
                    if line.find('percentiles of the PDF') >= 0:
                        percentile_mode = True
                        read_cols = [param_name+' P2.5', param_name+' P16', param_name+' P50', param_name+' P84', param_name+' P97.5']
                    else:
                        percentile_mode = False
                        param_name = 'PARAMETER ' + line
                        read_cols = [param_name+' BIN', param_name+' PDF']
            else:
                # this line is a non-commented data line
                line_split = line.split()
                if len(read_cols) > 0 and len(read_cols) == len(line_split):
                    read_vals = line_split
                    for read_col,read_val in itertools.zip_longest(read_cols,read_vals):
                        if read_col in fit_fit_data:
                            fit_fit_data[read_col].append(float(read_val))
                        else:
                            fit_fit_data[read_col] = [float(read_val)]
                        #print('fit_fit_data[\'%s\'] = %s' % (read_col, read_val)) # debug
                else:
                    continue
    # 
    return fit_fit_data



# 
# Function for converting energy units and scales
# 
def convert_energies_to_flux_densities(energies, energy_unit, wavelength_um = [], redshift = np.nan):
    SED_flux_mJy = []
    if re.search(r'\bLoA^-1\b', energy_unit):
        # if the energies are L_{\lambda} in units of L_{\odot} {\AA}^{-1}
        if len(wavelengths_um) == 0:
            print('******')
            print('Warning! convert_energies_to_flux_densities() requires \'wavelength_um\'! Conversion failed!')
            print('******')
        if np.isnan(redshift):
            print('******')
            print('Warning! convert_energies_to_flux_densities() requires \'redshift\'! Conversion failed!')
            print('******')
        vLv = energies * (wavelengths_um/1e4)
        lumdist_Mpc = cosmo.luminosity_distance(redshift).value # Mpc
        SED_flux_mJy = vLv / (4 * np.pi * lumdist_Mpc**2) * (1.+redshift) * 40.31970 / (2.99792458e5/(wavelength_um)) # 1 L_{\odot} Mpc^{-2} = 40.31970 mJy GHz
        #SED_Lv = vLv / (2.99792458e8/(wavelength_um/1e6)) # L_{\odot} Hz^{-1}
    elif re.search(r'\bJy\b', energy_unit):
        # if the energies are flux densities in units of Jy
        SED_flux_mJy = energies * 1e3
    # 
    if len(SED_flux_mJy) == 0:
        raise('Error!! convert_energies_to_flux_densities() failed!')
    # 
    return SED_flux_mJy



# 
# Function for splining flux densities with wavelenegths
# 
def spline_flux_densities(SED_flux_mJy, SED_wavelength_um, output_wavelength_um):
    _, index_unique = np.unique(SED_wavelength_um, return_index=True) # remove non-monochromatically-increasing data
    spliner = interpolate.UnivariateSpline(np.log10(SED_wavelength_um[index_unique]), np.log10(SED_flux_mJy[index_unique])) # spline in logarithm space
    output_flux_mJy = np.power(10, spliner(np.log10(output_wavelength_um ) ) )
    return output_flux_mJy



# 
# Subroutine of printing the usage
# 
def usage():
    print('Usage: ')
    print('  %s fit_name [-user-obs XXX.txt -filters XXX.txt]' % (os.path.basename(__file__) ) )



# 
# Main Program
# 
if __name__ == '__main__':
    
    obj_name = ''
    fit_names = []
    magphys_user_obs_file = 'magphys_input_fluxes.dat'
    magphys_filters_file = 'magphys_input_filters.dat'
    if 'USER_OBS' in os.environ:
        magphys_user_obs_file = os.environ('USER_OBS')
    if 'USER_FILTERS' in os.environ:
        magphys_filters_file = os.environ('USER_FILTERS')
    
    # Check user input
    if len(sys.argv) <= 1:
        usage()
        sys.exit()
    
    # Read user input
    i = 1
    while i < len(sys.argv):
        arg_str = sys.argv[i].lower().replace('--','-')
        if arg_str == '-help':
            usage()
            sys.exit()
        elif arg_str == '-user-obs':
            if i+1 < len(sys.argv):
                i = i+1
                magphys_user_obs_file = sys.argv[i]
        elif arg_str == '-filters':
            if i+1 < len(sys.argv):
                i = i+1
                magphys_user_obs_file = sys.argv[i]
        elif arg_str == '-source' or arg_str == '-obj' or arg_str == '-obj-name':
            if i+1 < len(sys.argv):
                i = i+1
                obj_name = sys.argv[i]
        else:
            fit_names.append(sys.argv[i])
        i = i+1
    
    # Check user input
    if len(fit_names) == 0:
        usage()
        sys.exit()
    
    # For each fit_name
    for fit_name in fit_names:
        if fit_name.endswith('.sed') or fit_name.endswith('.fit'):
            fit_name = fit_name[:-4]
        fit_sed_file = '%s.sed'%(fit_name)
        fit_fit_file = '%s.fit'%(fit_name)
        fit_res_file = '%s.res'%(fit_name) # residual flux
        fit_rf_file = '%s.rf'%(fit_name) # rest-frame flux
        fit_param_file = '%s.param'%(fit_name) # parameters
        if obj_name == '':
            obj_name = fit_name
        if os.path.isfile(fit_sed_file) and os.path.isfile(fit_fit_file):
            
            fit_sed_data = read_magphys_sed_file(fit_sed_file)
            #print(fit_sed_data)
            #-- Note that the best-fit SED flux in "*.sed" are in logarithm and have units of L_{\odot} AA^{-1}, not mJy! See magphys documentation.
            
            fit_fit_data = read_magphys_fit_file(fit_fit_file)
            #print(fit_fit_data)
            fit_fit_data_filters = [t for t in fit_fit_data if len(fit_fit_data[t])==3]
            #print(fit_fit_data_filters)
            #print(len(fit_fit_data_filters))
            #-- Note that the best-fit SED flux at each band in "*.fit" files have units of L_{\odot} Hz^{-1}, not mJy! See magphys documentation.
            
            user_obs_data = read_magphys_user_obs_file(magphys_user_obs_file)
            #-- Note that the observed flux at each band in "$USER_OBS" files have units of Jy, not mJy!
            
            filters_data = read_magphys_filters_file(magphys_filters_file)
            filter_names = filters_data[filters_data.colnames[0]]
            filter_lambdas = filters_data[filters_data.colnames[1]]
            filter_is_fit = filters_data[filters_data.colnames[3]]
            filter_obs_fluxes = []
            filter_obs_flux_errors = []
            #-- Note that the observed-frame wavelength at each band in "$USER_FILTERS" files (lambda_eff) have units of um!
            
            #-- Note that all the above SED-related wavelengths are observed frame!
            
            #-- Hereafter we assume that the user_obs file contains only one source, i.e., one data line
            
            # Check the consistency between USER_OBS and USER_FILTERS data
            if len(filters_data) != (len(user_obs_data.colnames) - 2) / 2:
                print('******')
                print('Error! Inconsistenct user_obs and user_filters files! Please make sure len(filters_data) == (len(user_obs_data.colnames) - 2) / 2!')
                print('******')
                sys.exit()
            # Rename USER_OBS columns because in some cases the column names have duplicates and could not be read by astropy.io.ascii well
            user_obs_data.rename_column(user_obs_data.colnames[0], 'obj_name')
            user_obs_data.rename_column(user_obs_data.colnames[1], 'redshift')
            for ifilter in range(len(filters_data)):
                # the first 2 columns in user_obs_data are obj_name and z, then the columns are flux and flux error at each filter_data band. 
                user_obs_data.rename_column(user_obs_data.colnames[ifilter*2+2], filter_names[ifilter])
                user_obs_data.rename_column(user_obs_data.colnames[ifilter*2+3], 'E_'+filter_names[ifilter])
                filter_obs_fluxes.append(user_obs_data[user_obs_data.colnames[ifilter*2+2]][0]) # assuming only one data line
                filter_obs_flux_errors.append(user_obs_data[user_obs_data.colnames[ifilter*2+3]][0]) # assuming only one data line
            #print(user_obs_data)
            #print(filters_data)
            filter_obs_fluxes = np.array(filter_obs_fluxes)
            filter_obs_flux_errors = np.array(filter_obs_flux_errors)
            
            
            
            # compute residual fluxes
            with open(fit_res_file, 'w') as ofp:
                ofp_fmt_width_of_obj_name = max( [ len(obj_name), len('obj_name') ] ) + 3
                ofp_fmt_width_of_filter_name = max( [ len(max(filter_names, key=len)), len('filter_name') ] ) + 2
                ofp_fmt = '# %%-%ds %%15s %%%ds %%10s %%15s %%15s %%15s %%15s %%15s\n' % (ofp_fmt_width_of_obj_name, ofp_fmt_width_of_filter_name )
                ofp.write(ofp_fmt % ('obj_name', 'redshift', 'filter_name', 'is_fit', 'wavelength_um', 'OBS_flux_mJy', 'E_OBS_flux_mJy', 'SED_flux_mJy', 'RES_flux_mJy') )
                for ifilter in range(len(filters_data)):
                    filter_name = filter_names[ifilter]
                    wavelength_um = filter_lambdas[ifilter] # um, obs-frame
                    is_fit = filter_is_fit[ifilter]
                    redshift = fit_fit_data['redshift'][0] # assuming only one data line
                    lumdist_Mpc = cosmo.luminosity_distance(redshift).value # Mpc
                    if filter_name in fit_fit_data:
                        if len(fit_fit_data[filter_name]) >= 3:
                            # convert L_{\odot} Hz^{-1} to mJy # <TODO> *.fit photometry data unit
                            luminosities = np.array(fit_fit_data[filter_name]) # L_{\odot} Hz^{-1}, should contain 3 values: observed luminosity, observed luminosity error, and SED best-fit luminosity. 
                            frequencies = 2.99792458e8/(wavelength_um/1e6) # Hz
                            vLv = luminosities * frequencies # L_{\odot}
                            fluxes_mJy = vLv / (4 * np.pi * lumdist_Mpc**2) * (1.+redshift) * 40.31970 / (2.99792458e5/(wavelength_um)) # 1 L_{\odot} Mpc^{-2} = 40.31970 mJy GHz
                            ofp_fmt = '  %%-%ds %%15.4f %%%ds %%10d %%15.6f %%15.6e %%15.6e %%15.6e %%15.6e\n' % (ofp_fmt_width_of_obj_name, ofp_fmt_width_of_filter_name )
                            ofp.write(ofp_fmt % (obj_name, redshift, filter_name, is_fit, wavelength_um, fluxes_mJy[0], fluxes_mJy[1], fluxes_mJy[2], fluxes_mJy[0]-fluxes_mJy[2]) ) # [0] is observed flux, [1] is observed flux error, [2] is SED best-fit flux.
                        else:
                            print('******')
                            print('Warning! "%s" filter %s data has less than 3 values?'%(fit_fit_file, filter_name) )
                            print('******')
                    else:
                        print('******')
                        print('Error! "%s" does not contain filter %s data?'%(fit_fit_file, filter_name) )
                        print('******')
                print('Written to "%s"' % (fit_res_file))
            
            
            # compute rest-frame fluxes with the best-fit SED from the SED fitting result "*.sed" file
            with open(fit_rf_file, 'w') as ofp:
                wavelengths = []
                energies = []
                # determine wavelengths and energies data
                if 'lg(lambda/A)' in fit_sed_data and 'Attenuated' in fit_sed_data:
                    wavelengths = np.power(10, np.array(fit_sed_data['lg(lambda/A)']) ) # original unit \AA, obs-frame, converted to um unit.
                    energies = np.power(10, np.array(fit_sed_data['Attenuated']) )
                    wavelength_unit = 'A' # determine wavelength unit
                    energy_unit = fit_sed_data['__energy_unit__'] # determine energy unit
                else:
                    # <TODO> if MAGPHYS *.sed file format changes, here we also need to apply some changes
                    print('******')
                    print('Error! "%s" does not contain "lg(lambda/A)" and "Attenuated" columns?'%(fit_fit_file, filter_name) )
                    print('******')
                # compute wavelength_um
                if len(wavelengths) > 0:
                    if wavelength_unit == 'A':
                        wavelength_um = wavelengths / 1e4 # um
                    else:
                        print('******')
                        print('Error! Could not determine wavelength unit!')
                        print('******')
                # compute SED_flux_mJy and spline the flux densities of rest-frame wavelengths
                if len(energies) > 0:
                    SED_flux_mJy = convert_energies_to_flux_densities(energies, energy_unit, wavelength_um, redshift)
                    RF_wavelength_um = np.array([850.0, 500.0, 350.0, 250.0])
                    RF_flux_mJy = spline_flux_densities(SED_flux_mJy, wavelength_um, RF_wavelength_um*(1.0+redshift))
                    #_, index_unique = np.unique(wavelength_um, return_index=True)
                    #spliner = interpolate.UnivariateSpline(np.log10(wavelength_um[index_unique]), np.log10(SED_flux_mJy[index_unique]))
                    #print('%0.6e' % np.power(10, spliner(np.log10(1297.32) ) ) ) # mJy, should be consistent with the values in the "$USER_OBS" file, but note that the latter file contains flux in units of Jy.
                    ##SED_Lv_spliner = interpolate.UnivariateSpline(np.log10(wavelength_um[index_unique]), np.log10(SED_Lv[index_unique]))
                    ###print('%0.6e' % np.power(10, SED_Lv_spliner(np.log10(1297.32) ) ) ) # Lsun Hz-1, should be consistent with the values in the "*.fit" file. 
                    for iRF in range(len(RF_wavelength_um)):
                        if iRF == 0:
                            ofp_fmt = '# %%-%ds %%15s %%18s %%18s %%16s\n' % (ofp_fmt_width_of_obj_name )
                            ofp.write(ofp_fmt % ('obj_name', 'redshift', 'RF_wavelength_um', 'OBS_wavelength_um', 'RF_flux_mJy' ) )
                        ofp_fmt = '  %%-%ds %%15.4f %%18.6f %%18.6f %%16.6e\n' % (ofp_fmt_width_of_obj_name )
                        ofp.write(ofp_fmt % (obj_name, redshift, RF_wavelength_um[iRF], RF_wavelength_um[iRF]*(1.0+redshift), RF_flux_mJy[iRF] ) )
                    # 
                    print('Written to "%s"' % (fit_rf_file))
            
            
            # write parameter result table
            with open(fit_param_file, 'w') as ofp:
                # we first compute SED fittin parameters
                #   nfit* is the number of fitted data points
                #   snr* is the S/N of fitted data points added in quadratic
                #   rchi2* is the reduced chi-square of fitted data points, rchi2 = \sum { (f_OBS - f_SED )^2 / e_OBS^2 } / ( Nfit - 1 )
                SED_fitting_parameters_name_dict = {}
                SED_fitting_parameters_name_dict['nfit'] = 'nfit'
                SED_fitting_parameters_name_dict['snr'] = 'snr'
                SED_fitting_parameters_name_dict['rchi2'] = 'rchi2'
                SED_fitting_parameters_value_dict = {}
                SED_fitting_parameters_value_dict['nfit'] = 0
                SED_fitting_parameters_value_dict['nfit_star'] = 0
                SED_fitting_parameters_value_dict['nfit_dust'] = 0
                SED_fitting_parameters_value_dict['snr'] = 0.0
                SED_fitting_parameters_value_dict['snr_star'] = 0.0
                SED_fitting_parameters_value_dict['snr_dust'] = 0.0
                SED_fitting_parameters_value_dict['rchi2'] = 0.0
                SED_fitting_parameters_value_dict['rchi2_star'] = 0.0
                SED_fitting_parameters_value_dict['rchi2_dust'] = 0.0
                for ifilter in range(len(filter_names)):
                    filter_name = filter_names[ifilter]
                    filter_lambda = filter_lambdas[ifilter]
                    is_fit = filter_is_fit[ifilter]
                    RF_wavelength_um = filter_lambda / (1.0+redshift) # convert from obs-frame to rest-frame
                    if fit_fit_data[filter_name][0]>0 and fit_fit_data[filter_name][1]>0 and is_fit == 1:
                        if RF_wavelength_um < 8.0:
                            if fit_fit_data[filter_name][0] > fit_fit_data[filter_name][1] * 3.0:
                                # S/N > 3 data points were fit
                                SED_fitting_parameters_value_dict['nfit'] = SED_fitting_parameters_value_dict['nfit'] + 1
                                SED_fitting_parameters_value_dict['nfit_star'] = SED_fitting_parameters_value_dict['nfit_star'] + 1
                                SED_fitting_parameters_value_dict['snr'] = SED_fitting_parameters_value_dict['snr'] + (fit_fit_data[filter_name][0] / fit_fit_data[filter_name][1])**2
                                SED_fitting_parameters_value_dict['snr_star'] = SED_fitting_parameters_value_dict['snr_star'] + (fit_fit_data[filter_name][0] / fit_fit_data[filter_name][1])**2
                                SED_fitting_parameters_value_dict['rchi2'] = SED_fitting_parameters_value_dict['rchi2'] + ((fit_fit_data[filter_name][0]-fit_fit_data[filter_name][2]) / fit_fit_data[filter_name][1])**2
                                SED_fitting_parameters_value_dict['rchi2_star'] = SED_fitting_parameters_value_dict['rchi2_star'] + ((fit_fit_data[filter_name][0]-fit_fit_data[filter_name][2]) / fit_fit_data[filter_name][1])**2
                        elif RF_wavelength_um < 2000.0:
                            if fit_fit_data[filter_name][0] > fit_fit_data[filter_name][1] * 2.0:
                                # S/N > 2 data points were fit
                                SED_fitting_parameters_value_dict['nfit'] = SED_fitting_parameters_value_dict['nfit'] + 1
                                SED_fitting_parameters_value_dict['nfit_dust'] = SED_fitting_parameters_value_dict['nfit_dust'] + 1
                                SED_fitting_parameters_value_dict['snr'] = SED_fitting_parameters_value_dict['snr'] + (fit_fit_data[filter_name][0] / fit_fit_data[filter_name][1])**2
                                SED_fitting_parameters_value_dict['snr_dust'] = SED_fitting_parameters_value_dict['snr_dust'] + (fit_fit_data[filter_name][0] / fit_fit_data[filter_name][1])**2
                                SED_fitting_parameters_value_dict['rchi2'] = SED_fitting_parameters_value_dict['rchi2'] + ((fit_fit_data[filter_name][0]-fit_fit_data[filter_name][2]) / fit_fit_data[filter_name][1])**2
                                SED_fitting_parameters_value_dict['rchi2_dust'] = SED_fitting_parameters_value_dict['rchi2_dust'] + ((fit_fit_data[filter_name][0]-fit_fit_data[filter_name][2]) / fit_fit_data[filter_name][1])**2
                SED_fitting_parameters_value_dict['snr'] = np.sqrt(SED_fitting_parameters_value_dict['snr'])
                SED_fitting_parameters_value_dict['snr_star'] = np.sqrt(SED_fitting_parameters_value_dict['snr_star'])
                SED_fitting_parameters_value_dict['snr_dust'] = np.sqrt(SED_fitting_parameters_value_dict['snr_dust'])
                if float(SED_fitting_parameters_value_dict['nfit'] - 1) > 0:
                    SED_fitting_parameters_value_dict['rchi2'] = (SED_fitting_parameters_value_dict['rchi2']) / float(SED_fitting_parameters_value_dict['nfit'] - 1)
                else:
                    SED_fitting_parameters_value_dict['rchi2'] = np.nan
                if float(SED_fitting_parameters_value_dict['nfit_star'] - 1) > 0:
                    SED_fitting_parameters_value_dict['rchi2_star'] = (SED_fitting_parameters_value_dict['rchi2_star']) / float(SED_fitting_parameters_value_dict['nfit_star'] - 1)
                else:
                    SED_fitting_parameters_value_dict['rchi2_star'] = np.nan
                if float(SED_fitting_parameters_value_dict['nfit_dust'] - 1) > 0:
                    SED_fitting_parameters_value_dict['rchi2_dust'] = (SED_fitting_parameters_value_dict['rchi2_dust']) / float(SED_fitting_parameters_value_dict['nfit_dust'] - 1)
                else:
                    SED_fitting_parameters_value_dict['rchi2_dust'] = np.nan
                
                # then we prepare to read SED fitting derived physical parameters from 'fit_fit_data'
                SED_physical_parameters_name_dict = {}
                SED_physical_parameters_name_dict['logMstar'] = 'M(stars)'
                SED_physical_parameters_name_dict['logMdust'] = 'M(dust)'
                SED_physical_parameters_name_dict['logLdust'] = 'Ldust'
                SED_physical_parameters_name_dict['logSFR'] = 'SFR_0.1Gyr'
                SED_physical_parameters_name_dict['logsSFR'] = 'sSFR_0.1Gyr'
                SED_physical_parameters_name_dict['logAge'] = 'age_M (mass-weighted age)'
                SED_physical_parameters_name_dict['A_V'] = 'A_V'
                SED_physical_parameters_name_dict['Tdust'] = 'Tdust'
                fit_fit_data['M(stars)'] = np.log10(fit_fit_data['M*'])
                fit_fit_data['M(dust)'] = np.log10(fit_fit_data['Mdust'])
                fit_fit_data['Ldust'] = np.log10(fit_fit_data['Ldust'])
                fit_fit_data['SFR_0.1Gyr'] = np.log10(fit_fit_data['SFR'])
                fit_fit_data['sSFR_0.1Gyr'] = np.log10(fit_fit_data['sSFR'])
                fit_fit_data['age_M (mass-weighted age)'] = fit_fit_data['age_M']
                
                # -- write obj_name and reshift column names into the header line
                ofp_fmt = '# %%-%ds %%9s' % (ofp_fmt_width_of_obj_name )
                ofp.write(ofp_fmt % ('obj_name', 'redshift' ) )
                # -- write nfit nfit_star nfit_dust chisq chisq_star chisq_dust column names into the header line
                for par_name in SED_fitting_parameters_name_dict:
                    ofp_fmt = ' %%%ds' % (max( [ len(par_name)+5, 9 ] ) ) # optimize the width of column
                    ofp.write(ofp_fmt % (par_name ) )
                    ofp.write(ofp_fmt % (par_name+'_star' ) )
                    ofp.write(ofp_fmt % (par_name+'_dust' ) )
                # -- write paramter column headers in the header line
                for par_name in SED_physical_parameters_name_dict:
                    ofp_fmt = ' %%%ds' % (max( [ len(par_name)+4, 7 ] ) ) # optimize the width of column
                    ofp.write(ofp_fmt % (par_name ) )
                    ofp.write(ofp_fmt % (par_name+'_MED' ) )
                    ofp.write(ofp_fmt % (par_name+'_H68' ) )
                    ofp.write(ofp_fmt % (par_name+'_L68' ) )
                # -- finished writing the header line
                ofp.write('\n')
                # -/ write obj_name and reshift into the data line(s)
                ofp_fmt = '  %%-%ds %%9.4f' % (ofp_fmt_width_of_obj_name )
                ofp.write(ofp_fmt % (obj_name, redshift ) )
                # -- write nfit nfit_star nfit_dust chisq chisq_star chisq_dust values into the data line(s)
                for par_name in SED_fitting_parameters_name_dict:
                    ofp_fmt = ' %%%dg' % (max( [ len(par_name)+5, 9 ] ) ) # optimize the width of column
                    ofp.write(ofp_fmt % (SED_fitting_parameters_value_dict[par_name] ) )
                    ofp.write(ofp_fmt % (SED_fitting_parameters_value_dict[par_name+'_star'] ) )
                    ofp.write(ofp_fmt % (SED_fitting_parameters_value_dict[par_name+'_dust'] ) )
                # -/ write parameter fit values into the data line(s)
                for par_name in SED_physical_parameters_name_dict:
                    ofp_fmt = ' %%%d.3f' % (max( [ len(par_name)+4, 7 ] ) ) # optimize the width of column
                    ofp.write(ofp_fmt % (fit_fit_data[SED_physical_parameters_name_dict[par_name]][0] ) ) # best-fit value from the "*.fit" file
                    ofp.write(ofp_fmt % (fit_fit_data['PARAMETER ' + SED_physical_parameters_name_dict[par_name] + ' P50'][0] ) )
                    ofp.write(ofp_fmt % (fit_fit_data['PARAMETER ' + SED_physical_parameters_name_dict[par_name] + ' P84'][0] ) )
                    ofp.write(ofp_fmt % (fit_fit_data['PARAMETER ' + SED_physical_parameters_name_dict[par_name] + ' P16'][0] ) )
                ofp.write('\n')
                # -- finished writing the data line(s)
                # 
                print('Written to "%s"' % (fit_param_file))
            
            
        else:
            print('******')
            print('Warning! %s or %s does not exist!' % (fit_sed_file, fit_fit_file))
            print('******')
















 
 
 