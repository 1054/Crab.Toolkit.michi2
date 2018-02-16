#!/usr/bin/env python
# 

# 
#   extracted_magnitude_Arp193.csv
#   extracted_magnitude_Arp220.csv
#   extracted_magnitude_IRASF17207-0014.csv
#   extracted_magnitude_NGC1614.csv
#   extracted_magnitude_NGC2623.csv
#   extracted_magnitude_NGC6240.csv
#   extracted_magnitude_UGC05101.csv
# 


import os, sys, glob, math, astropy
import astropy.io.ascii as asciitable

photometry_files_in_mag = glob.glob('datatable_photometry/Brown2014/extracted_magnitude_*.csv')
#photometry_files_in_mag = glob.glob('datatable_photometry/Brown2014/extracted_magnitude_Arp220_2.csv') #<DEBUG>#

wavelength_dict = {}
wavelength_dict['GALEX_FUV']     = 0.1531        # FUV           # 0.1531
wavelength_dict['Swift_UVW2']    = 0.2026        # UVW2          # 0.2026
wavelength_dict['Swift_UVM2']    = 0.2238        # UVM2          # 0.2238
wavelength_dict['GALEX_NUV']     = 0.2286        # NUV           # 0.2286
wavelength_dict['Swift_UVM1']    = 0.2598        # UVW1          # 0.2598
wavelength_dict['Swift_U']       = 0.3459        # U             # 0.3459
wavelength_dict['SDSS_u']        = 0.3551        # u             # 0.3551
wavelength_dict['SDSS_g']        = 0.4681        # g             # 0.4681
wavelength_dict['Swift_V']       = 0.5419        # V             # 0.5419
wavelength_dict['SDSS_r']        = 0.6165        # r             # 0.6165
wavelength_dict['SDSS_i']        = 0.7480        # i             # 0.7480
wavelength_dict['SDSS_z']        = 0.8931        # z             # 0.8931
wavelength_dict['2MASS_J']       = 1.232         # J             # 1.232
wavelength_dict['2MASS_H']       = 1.644         # H             # 1.644
wavelength_dict['2MASS_Ks']      = 2.159         # KS            # 2.159
wavelength_dict['WISE_W1']       = 3.357         # W1            # 3.357
wavelength_dict['IRAC_I1']       = 3.544         # [3.6]         # 3.544
wavelength_dict['IRAC_I2']       = 4.487         # [4.5]         # 4.487
wavelength_dict['WISE_W2']       = 4.606         # W2            # 4.606
wavelength_dict['IRAC_I3']       = 5.710         # [5.8]         # 5.710
wavelength_dict['IRAC_I4']       = 7.841         # [8.0]         # 7.841
wavelength_dict['WISE_W3']       = 11.81         # W3            # 11.81
wavelength_dict['IRS_PB']        = 15.80         # PUI blue      # 15.80
wavelength_dict['WISE_W4']       = 22.14         # W4            # 22.14
wavelength_dict['WISE_W42']      = 22.14         # W4′           # 22.14
wavelength_dict['IRS_PR']        = 22.32         # PUI red       # 22.32
wavelength_dict['MIPS_M1']       = 23.51         # [24]          # 23.51

output_dir = 'SED_fitting_michi2'
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

for photometry_file_in_mag in photometry_files_in_mag:
    # 
    wavelength_array = []
    flux_array = []
    fluxerr_array = []
    # 
    photometry_data = asciitable.read(photometry_file_in_mag, header_start=0, data_start=1, delimiter=',', guess=False)
    galaxy_name = os.path.basename(photometry_file_in_mag).replace('extracted_magnitude_','').replace('.csv','')
    if not os.path.isdir(output_dir + os.sep + galaxy_name):
        os.makedirs(output_dir + os.sep + galaxy_name)
    output_file = output_dir + os.sep + galaxy_name + os.sep + 'extracted_flux_Brown2014.txt'
    print(photometry_data)
    print('')
    for x in photometry_data.colnames:
        if x!='Name' and not x.endswith('_err'):
            #print(x)
            mag_value = float(photometry_data.field(x)[0])
            magerr_value = float(photometry_data.field(x+'_err')[0])
            flux_value = math.pow(10, mag_value / (-2.5)) * 3631 * 1e3 # mJy
            fluxerr_value = magerr_value * flux_value
            wavelength_array.append(wavelength_dict[x])
            flux_array.append(flux_value)
            fluxerr_array.append(fluxerr_value)
    asciitable.write([wavelength_array, flux_array, fluxerr_array], sys.stdout, Writer=asciitable.FixedWidthTwoLine, delimiter='|', bookend=True, 
                     names=['wavelength_um','flux_mJy','e_flux_mJy'], 
                     formats={'wavelength_um':'%20.10f','flux_mJy':'%20.10f','e_flux_mJy':'%20.10f'})
    asciitable.write([wavelength_array, flux_array, fluxerr_array], output_file, Writer=asciitable.FixedWidthTwoLine, delimiter=' ', bookend=True, 
                     names=['wavelength_um','flux_mJy','e_flux_mJy'], 
                     formats={'wavelength_um':'%20.10f','flux_mJy':'%20.10f','e_flux_mJy':'%20.10f'}, overwrite=True)
    print('')
    print('Output to "%s"!'%(output_file))
    print('')
    #break



