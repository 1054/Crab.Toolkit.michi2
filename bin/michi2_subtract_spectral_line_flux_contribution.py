#!/usr/bin/env python
# 
# This code will subtract spectral line flux contribution to the SED photometry flux densities. 
# Inputs: redshift, 
# 
# 

import os, sys, re, json
import numpy as np
import astropy
from astropy.table import Table, Column
from copy import copy
import shutil


sys.path.append('/Users/dzliu/Cloud/Github/Crab.Toolkit.Python/lib/crab/crabpdbi')
from CrabPdBI import convert_Wavelength_um_to_Frequency_GHz, \
                    calc_radio_line_frequencies, \
                    find_radio_lines_in_frequency_range, \
                    calc_radio_line_flux_from_IR_luminosity




#####################################
#               USAGE               #
#####################################

def usage():
    print('Usage: ')
    print('    michi2_subtract_spectral_line_flux_contribution.py datatable_photometry.txt \\')
    print('                                                       -output datatable_photometry_output.txt \\')
    print('                                                       -redshift 3.0 \\')
    print('                                                       -IR-luminosity 5e12 \\')
    print('                                                       -freq-support extracted_freq_support.json \\')
    print('                                                       [-line-names CO CI NII CII] \\')
    print('                                                       [-wavelength-column 1] \\')
    print('                                                       [-flux-value-column 2] \\')
    print('                                                       [-flux-error-column 3] \\')
    print('                                                       [-flux-unit-column 4] \\')
    print('                                                       [-filter-name-column 5] \\')
    print('                                                       [-known-line-name-and-flux \"CO(2-1)\" \"3.0 Jy km/s\"] \\')
    print('                                                       [-known-line-name-and-flux \"CO(3-2)\" \"5.0 Jy km/s\"] \\')
    print('                                                       ')
    print('Notes:')
    print('    The first input must be a text format photometry data table, e.g., named "datatable_photometry.txt".')
    print('    The photometry data table "datatable_photometry.txt" should have 5 columns: wavelength (um), flux, flux error, unit, and filter name.')
    print('    The input -filter-freq-support should be a json format file containing a Python dictionary, with keys being filter names and values being frequency lists.')
    print('    TODO')



####################################
#               MAIN               #
####################################

if __name__ == '__main__':
    
    # Read user input
    input_table_file = ''
    input_redshift = np.nan
    input_freq_support_json_files = []
    input_freq_support_filter_names = []
    input_freq_support_freq_lists = []
    input_IR_luminosity = np.nan
    input_line_names = ['CO','CI','NII','CII']
    input_filter_names = []
    input_filter_regexes = []
    input_wavelength_column = 1
    input_flux_value_column = 2
    input_flux_error_column = 3
    input_flux_unit_column = 4
    input_filter_name_column = 5
    known_line_names = []
    known_line_fluxes = []
    output_file_name = ''
    speed_of_light_kms = 2.99792458e5
    iarg = 1
    tstr = ''
    tmode = ''
    isopt = False
    while iarg < len(sys.argv):
        tstr = re.sub(r'^[-]+', r'-', sys.argv[iarg].lower()) # lower case current input argument string
        isopt = False # whether current input argument is an option (starting with "-")
        if tstr == '-redshift' or tstr == '-z':
            tmode = 'redshift'
            isopt = True
        elif tstr == '-IR-luminosity'.lower():
            tmode = 'IR_luminosity'
            isopt = True
        elif tstr == '-freq-support':
            tmode = 'freq_support'
            isopt = True
        elif tstr == '-line-names':
            tmode = 'line_names'
            isopt = True
        elif tstr == '-known-line-name-and-flux':
            tmode = 'known_lines'
            isopt = True
        elif tstr == '-wavelength-column':
            tmode = 'wavelength_column'
            isopt = True
        elif tstr == '-filter-name-column':
            tmode = 'filter_name_column'
            isopt = True
        elif tstr == '-output-file-name' or tstr == '-output-file' or tstr == '-output' or tstr == '-out':
            tmode = 'output_file_name'
            isopt = True
        elif tstr.startswith('-'):
            raise ValueError('The input argument "%s" is not allowed!'%(sys.argv[iarg]))
            sys.exit()
            #<TODO># what if some value starts with "-"?
        else:
            if tmode == '':
                input_table_file = sys.argv[iarg]
        # 
        #print(sys.argv[iarg])
        #print('isopt', isopt)
        #print('tmode', tmode)
        # 
        if isopt == False:
            if tmode == 'redshift':
                input_redshift = float(sys.argv[iarg])
            elif tmode == 'IR_luminosity':
                input_IR_luminosity = float(sys.argv[iarg])
            elif tmode == 'freq_support':
                if sys.argv[iarg].endswith('.json'):
                    input_freq_support_json_files.append(sys.argv[iarg])
                    with open(sys.argv[iarg], 'r') as fp:
                        freq_dict = json.load(fp)
                        for key in freq_dict:
                            input_freq_support_filter_names.append(key)
                            input_freq_support_freq_lists.append(freq_dict[key])
                else:
                    try:
                        if iarg + 1 < len(sys.argv):
                            input_freq_support_filter_names.append(sys.argv[iarg])
                            input_freq_support_freq_lists.append(eval('['+sys.argv[iarg+1]+']'))
                            iarg = iarg + 1
                        else:
                            raise ValueError('The input freq_support "%s" is neither a json file nor a filter_name freq_list pair!'%(sys.argv[i]))
                            sys.exit()
                    except:
                        raise ValueError('The input freq_support "%s" is neither a json file nor a filter_name freq_list pair!'%(sys.argv[i]))
                        sys.exit()
            elif tmode == 'line_names':
                input_line_names.append(sys.argv[iarg])
            elif tmode == 'filter_names':
                input_filter_names.append(sys.argv[iarg])
                input_filter_regexes.append(re.compile(sys.argv[iarg]))
            elif tmode == 'wavelength_column':
                input_wavelength_column = int(sys.argv[iarg])
            elif tmode == 'filter_name_column':
                input_filter_name_column = int(sys.argv[iarg])
            elif tmode == 'known_lines':
                if iarg + 1 < len(sys.argv):
                    known_line_names.append(sys.argv[iarg])
                    known_line_fluxes.append(sys.argv[iarg+1])
                    iarg = iarg + 1
            elif tmode == 'output_file_name':
                output_file_name = sys.argv[iarg]
        # 
        iarg = iarg + 1
    
    # Check user input
    if input_table_file == '' or np.isnan(input_redshift) or len(input_freq_support_filter_names) == 0 or len(input_freq_support_freq_lists) == 0 or np.isnan(input_IR_luminosity) or output_file_name == '':
        usage()
        sys.exit()
    
    # Store variables
    z = input_redshift
    IR_luminosity = input_IR_luminosity
    
    # Read data table
    tb = Table.read(input_table_file, format='ascii')
    tbout = tb.copy()
    
    # Loop each photometry
    wavelengths = tb.field(tb.colnames[input_wavelength_column-1]).data # um
    fluxvalues = tb.field(tb.colnames[input_flux_value_column-1]).data # mJy
    fluxerrors = tb.field(tb.colnames[input_flux_error_column-1]).data # mJy
    fluxunits = tb.field(tb.colnames[input_flux_unit_column-1]).data # str
    filternames = tb.field(tb.colnames[input_filter_name_column-1]).data # str
    frequencies = speed_of_light_kms / wavelengths # GHz
    subtractedlinename = []
    subtractedlineflux_Jykms = []
    subtractedlineflux_mJy = []
    countlines = 0
    for i in range(len(filternames)):
        # 
        subtracted_line_name = []
        subtracted_line_flux_Jykms = []
        subtracted_line_flux_mJy = []
        new_flux_value = np.nan
        new_flux_error = np.nan
        # 
        if filternames[i] in input_freq_support_filter_names:
            print('-----------------------------------------------------------------')
            filter_name = filternames[i]
            filter_flux = fluxvalues[i]
            filter_ferr = fluxerrors[i]
            original_flux = filter_flux
            original_ferr = filter_ferr
            print('filter_name:', "filter_name")
            print('filter_flux:', filter_flux)
            freq_list = input_freq_support_freq_lists[input_freq_support_filter_names.index(filternames[i])]
            filter_freqcenter = np.mean(np.array(freq_list))
            filter_bandwidth = np.sum(np.array(freq_list)[1::2]-np.array(freq_list)[0::2])
            print('freq_list:', freq_list)
            print('filter_freqcenter:', filter_freqcenter)
            print('filter_bandwidth:', filter_bandwidth)
            line_freqs, line_names = find_radio_lines_in_frequency_range(freq_list, Redshift=z, include_faint_lines = False)
            # 
            # -- 20181203 bug fixed: there are some duplicates in line_names if the spectral_setups have overlaps
            j = len(line_names)-1
            while j > 0:
                if line_names[j] in line_names[0:j]:
                    del line_names[j]
                    del line_freqs[j]
                j -= 1
            # 
            # Loop line names within frequency support
            for j in range(len(line_names)):
                line_name = line_names[j]
                line_freq = line_freqs[j]
                # 
                # predict line flux
                line_flux_prediction = calc_radio_line_flux_from_IR_luminosity(line_name, IR_luminosity, z, verbose=True) # TODO: IR_color = IR_color
                # 
                # 20180830 -- get real CO observation
                #for k in range(len(radio_line_database_table)):
                #    database_line_freq, database_line_name = calc_radio_line_frequencies(radio_line_database_table['LineName'][k], set_output_line_names = True) # recognize the line name 
                #    if radio_line_database_table['ID'][k] == ID_prior and database_line_name == line_name:
                #        output_line_prediction_vs_observation_fp.write('  %-15d %18s %25s %18.4f %18.4f %18.4f %18s\n' % (i, ID_prior, line_name, line_flux_prediction, radio_line_database_table['LineFlux'][k], radio_line_database_table['LineFluxError'][k], radio_line_database_table['Source'][k]) )
                #        print('line_flux_density_in_databse = %0.4f +- %0.4f (vs. prediction %0.4f) [Jy km s-1] (ID %d) (Line Name "%s") (Source Name "%s")' % (radio_line_database_table['LineFlux'][k], radio_line_database_table['LineFluxError'][k], line_flux_prediction, ID_prior, line_name, radio_line_database_table['Source'][k] ) )
                # 
                # line_flux_contribution_as_a_continuum
                line_flux_contribution_as_a_continuum = line_flux_prediction / (filter_bandwidth/filter_freqcenter*speed_of_light_kms) * 1e3 # mJy, line flux distributed over 8 GHz at filter_freqcenter GHz
                print('line_flux_contribution_as_a_continuum = %0.4f [mJy]' % (line_flux_contribution_as_a_continuum))
                line_flux_contribution_to_the_continuum = line_flux_contribution_as_a_continuum / filter_flux
                print('line_flux_contribution_to_the_measured_flux = %0.4f%%' % (line_flux_contribution_to_the_continuum * 100))
                # 
                # 20180828: we need iterate a little bit, because IR_luminosity is based on the measured flux by assuming no line contamination
                # if there has a line contamination, say 
                #   S_line_as_a_continuum =    a  * S_measured_flux
                #   S_pure_dust_continuum = (1-a) * S_measured_flux
                # we also naively assume that IR_luminosity is linearly proportional to the pure dust continuum flux
                #   L_measured_IR_luminosity = b * S_measured_flux
                #       L_pure_IR_luminosity = b * S_pure_dust_continuum = (1-a) * L_measured_IR_luminosity
                # and we have the correlation between IR_luminosity and S_line_as_a_continuum
                #   S_line_as_a_continuum = f(IR_luminosity, z)
                # then we have
                #   f(L_measured_IR_luminosity, z) = a0 * S_measured_flux
                # and
                #   f(L_pure_IR_luminosity, z) = a * S_measured_flux
                #   i.e., f((1-a) * L_measured_IR_luminosity, z) = a * S_measured_flux
                # 
                print('Iterating to compute the line flux contribution fraction ...')
                a = 0.0
                while True:
                    line_flux_prediction = calc_radio_line_flux_from_IR_luminosity(line_name, (1-a) * IR_luminosity, z, verbose=False) # Jy km/s
                    line_flux_contribution_as_a_continuum = line_flux_prediction / (filter_bandwidth/filter_freqcenter*speed_of_light_kms) * 1e3 # mJy, line flux distributed over 8 GHz at filter_freqcenter GHz
                    print(line_flux_contribution_as_a_continuum, a * filter_flux, a)
                    if abs(line_flux_contribution_as_a_continuum - a * filter_flux) < 1e-6 * filter_flux:
                        break
                    a = line_flux_contribution_as_a_continuum / filter_flux
                    if a >= 1:
                        a = 1
                        break
                line_flux_contribution_to_the_continuum = a
                print('line_flux_contribution_to_the_measured_flux = %0.4f%% (after iteration)' % (line_flux_contribution_to_the_continuum * 100))
                # 
                # decontamination
                line_flux_contribution_as_a_continuum = line_flux_contribution_to_the_continuum * filter_flux
                filter_flux = (1-line_flux_contribution_to_the_continuum) * filter_flux # if there are multiple lines, this can iteratively remove all line contaminations
                filter_ferr = original_ferr # keep original flux error rather than S/N
                #filter_ferr = np.sqrt((original_ferr/original_flux)**2 ) * filter_flux # keep original S/N
                #filter_ferr = np.sqrt((original_ferr/original_flux)**2 + len(line_names) * 0.2**2) * filter_flux # assuming each line prediction has an error of 0.2 dex, then the error propagation is as this line of code.
                # 
                # print message
                print('"line_name":"%s", z:%0.6f, "IR_luminosity_for_prediction":%0.6e, "line_flux_prediction":%0.12f, '%(line_name, z, (1-line_flux_contribution_to_the_continuum) * IR_luminosity, line_flux_prediction))
                print('"line_flux_as_a_continuum":%0.12f, "total_filter_flux_measured":%0.12f, "filter_flux_after_line_subtraction":%0.12f, "ratio_of_line_contribution":%0.12f, '%(line_flux_contribution_as_a_continuum, original_flux, filter_flux, line_flux_contribution_to_the_continuum))
                subtracted_line_name.append(line_name)
                subtracted_line_flux_Jykms.append('%0.4e'%(line_flux_prediction))
                subtracted_line_flux_mJy.append('%0.4e'%(line_flux_contribution_as_a_continuum))
                
                countlines += 1
        # 
        if len(subtracted_line_name) > 0:
            subtractedlinename.append('"'+';'.join(subtracted_line_name)+'"')
            subtractedlineflux_Jykms.append('"'+';'.join(subtracted_line_flux_Jykms)+'"')
            subtractedlineflux_mJy.append('"'+';'.join(subtracted_line_flux_mJy)+'"')
            fluxvalues[i] = filter_flux
            fluxerrors[i] = filter_ferr
        else:
            subtractedlinename.append('""')
            subtractedlineflux_Jykms.append('""')
            subtractedlineflux_mJy.append('""')
    # 
    # 
    # save to output table
    tbout[tb.colnames[input_wavelength_column-1]] = wavelengths # um
    tbout[tb.colnames[input_flux_value_column-1]] = fluxvalues # mJy
    tbout[tb.colnames[input_flux_error_column-1]] = fluxerrors # mJy
    tbout[tb.colnames[input_flux_unit_column-1]] = fluxunits # str
    tbout[tb.colnames[input_filter_name_column-1]] = filternames # str
    tbout.add_column(Column(np.array(subtractedlinename)), name='subt_line_name')
    tbout.add_column(Column(np.array(subtractedlineflux_Jykms)), name='subt_line_flux_Jykms')
    tbout.add_column(Column(np.array(subtractedlineflux_mJy)), name='subt_line_flux_mJy')
    if os.path.isfile(output_file_name):
        print('Found existing "%s"! Backup as "%s.backup"!'%(output_file_name, output_file_name))
        shutil.move(output_file_name, output_file_name+'.backup')
    tbout.write(output_file_name, format='ascii.fixed_width', delimiter='  ', bookend=True, overwrite=False)
    with open(output_file_name, 'r+') as fp:
        fp.seek(0)
        fp.write('#')
    print('Output to "%s"!'%(output_file_name))
    print('')
    # 







