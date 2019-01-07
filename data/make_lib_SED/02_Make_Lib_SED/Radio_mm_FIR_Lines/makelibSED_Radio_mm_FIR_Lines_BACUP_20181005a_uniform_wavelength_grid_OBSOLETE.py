#!/usr/bin/env python
# 
# run this after running "a_dzliu_code_output_datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_SED_fluxes.py"
# 

from __future__ import print_function
import pkg_resources
pkg_resources.require('astropy>=2.0') # since 2.0 astropy.modeling.blackbody
import os, sys, re, json, time, astropy
import numpy as np
import astropy.io.ascii as asciitable
from astropy.table import Table, Column, hstack
#from astropy.convolution import convolve_fft
from numpy.fft import fft, ifft, fftfreq, fftshift, ifftshift
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from copy import copy
from pprint import pprint
from datetime import datetime
#from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)
from astropy import units as u
from astropy import constants as const
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
import scipy
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = '12' # https://matplotlib.org/users/customizing.html
mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
#mpl.rcParams['grid.color'] = 'b0b0b0'
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.linewidth'] = 0.25
mpl.rcParams['grid.alpha'] = 0.8




sys.path.append('/Users/dzliu/Cloud/Github/Crab.Toolkit.Python/lib/crab/crabpdbi/')
from CrabPdBI import ( find_radio_lines_in_frequency_range, 
                       calc_radio_line_flux_from_IR_luminosity, 
                       convert_flux2lprm, 
                       convert_flux2lsun, 
                     )





def get_SED_Radio_mm_FIR_lines(L_FIR = 1e12, linewidth = 500.0, starburst = 0.0, z = 0.0, verbose = False):
    # 
    # setup wavelength grid
    rest_wavelength_um_log10 = np.arange(np.log10(40.0), np.log10(1e5)+0.05, 0.0001) # 3um to 1e6um, grid 0.001 in log
    rest_wavelength_um = np.power(10.0, rest_wavelength_um_log10)
    nu = 2.99792458e5/rest_wavelength_um # GHz
    print('get_SED_Radio_mm_FIR_lines() grid number %d'%(len(rest_wavelength_um_log10)))
    # 
    # get all lines
    output_line_freqs, output_line_names = find_radio_lines_in_frequency_range([np.min(nu), np.max(nu)], Redshift = 0.0, set_output_line_names = True, include_faint_lines = False)
    rest_Sv = nu * 0.0*u.mJy
    for i in range(len(output_line_freqs)):
        IR_color = 0.6 + starburst #<TODO># IR_color
        IR_luminosity = L_FIR
        line_name = output_line_names[i]
        rest_freq = output_line_freqs[i]
        dL = cosmo.luminosity_distance(z).to(u.Mpc).value
        line_flux = calc_radio_line_flux_from_IR_luminosity(line_name, IR_luminosity, z, starburstiness = starburst, IR_color = IR_color, verbose = False)
        line_lprm = convert_flux2lprm(line_flux, rest_freq, z)
        # 
        # now this line_flux is computed at dL and z now. 
        # we need to remove these dependencies to construct the model. 
        # 
        # since line_lprm scales with (line_flux * dL**2 / (1.+z)) and is intrinsic, 
        # if in SED fitting, we use this model and fit a normalization coefficient 'aCOE'
        # which makes:
        #     line_flux_obs = aCOE * (line_flux)
        # then 
        #     (line_lprm_obs / dL_obs**2 * (1.+z_obs)) = aCOE * (line_lprm / dL**2 * (1.+z))
        # then 
        #      line_lprm_obs = aCOE * dL_obs**2 / (1.+z_obs) * line_lprm / dL**2 * (1.+z)
        # 
        # we renormalize line_flux to dL = 1 Mpc and (1.+z), 
        #     line_flux_SdLz = line_flux / dL**2 * (1.+z)
        # consequently
        #     line_lprm_SdLz = line_lprm / dL**2 * (1.+z)
        # then
        #     line_lprm_obs = aCOE * dL_obs**2 / (1.+z_obs) * line_lprm_SdLz
        # so
        #     IR_luminosity_obs = aCOE * dL_obs**2 / (1.+z_obs) * IR_luminosity_SdLz
        # 
        # In SED fitting, 
        # once we fit a normalization coefficient 'aCOE' to this SED model, 
        # then we can get 
        #     IR_luminosity_obs = aCOE * dL_obs**2 / (1.+z_obs) # unit is [1e12 Lsun]
        # 
        #####line_flux_SdLz = line_flux / dL**2 * (1.+z) # Jy km/s # must make sure cosmo here is the same as used in CrabPdBI.py
        #####IR_luminosity_SdLz = IR_luminosity / dL**2 * (1.+z)
        # 
        obs_freq = rest_freq / (1.+z)
        line_vel = (nu - rest_freq) / rest_freq * 2.99792458e5 # km/s
        line_sig = linewidth / (2.0*np.sqrt(2.0*np.log(2.0))) # km/s, from Gaussian FWHM to Gaussian sigma
        line_spec = 1.0/(line_sig*np.sqrt(2.0*np.pi)) * np.exp(-(line_vel)**2/(2.0*line_sig**2))
        line_spec = line_spec * line_flux*1e3 # mJy, area of the spec is line_flux*1e3 [mJy km/s].
        # 
        rest_Sv = rest_Sv + line_spec*u.mJy
        # 
        #print('Line %s with flux %s [Jy km/s] at obs freq %s [GHz] rest freq %s [GHz]'%(line_name, line_flux, obs_freq, rest_freq))
        print('Line %s with lprm %e [K km s-1 pc2] flux %s [Jy km/s] peak %s [mJy] at obs freq %s [GHz] rest freq %s [GHz]'%(line_name, line_lprm, line_flux, np.max(line_spec), obs_freq, rest_freq))
    # 
    # 
    if z > 0.0:
        obs_wavelength_um = rest_wavelength_um * (1.+z)
        obs_Sv = rest_Sv * (1.+z) # * dL**2 / (1.+z) # rest_Sv is normalized to L_FIR = 1e12, so here we multiply the real L_FIR.
    # 
    SED_dict = {}
    SED_dict['lambda_um'] = rest_wavelength_um
    SED_dict['Snu_mJy'] = rest_Sv.to(u.mJy).value # convert Jy to mJy
    SED_dict['Snu_Jy'] = rest_Sv.to(u.Jy).value
    SED_dict['lambda_rest_um'] = rest_wavelength_um
    SED_dict['Snu_rest_mJy'] = rest_Sv.to(u.mJy).value # convert Jy to mJy
    SED_dict['Snu_rest_Jy'] = rest_Sv.to(u.Jy).value
    if z > 0.0:
        SED_dict['lambda_obs_um'] = obs_wavelength_um
        SED_dict['Snu_obs_mJy'] = obs_Sv.to(u.mJy).value # convert Jy to mJy
        SED_dict['Snu_obs_Jy'] = obs_Sv.to(u.Jy).value
    return SED_dict



def get_Herschel_BandPass():
    BandPass_dict = {}
    BandPass_dict['PACS_70'] = {}
    BandPass_dict['PACS_100'] = {}
    BandPass_dict['PACS_160'] = {}
    BandPass_dict['SPIRE_250'] = {}
    BandPass_dict['SPIRE_350'] = {}
    BandPass_dict['SPIRE_500'] = {}
    for key in BandPass_dict:
        BandPass_dict[key]['lambda_um'] = []
        BandPass_dict[key]['response'] = []
        with open(os.path.dirname(os.path.abspath(__file__))+os.sep+'filters/Herschel/%s.dat'%(key), 'r') as fp:
            line = ''
            for line in fp:
                if line.startswith('#'): continue
                t = line.split()
                if len(t) < 2: continue
                BandPass_dict[key]['lambda_um'].append(float(t[0])/1e4) # convert AA to um
                BandPass_dict[key]['response'].append(float(t[1]))
    # 
    return BandPass_dict



def plot_Herschel_BandPass(ax):
    BandPass_dict = get_Herschel_BandPass()
    # 
    for key in BandPass_dict:
        x_filter_array = np.array(BandPass_dict[key]['lambda_um'])
        y_filter_array = np.array(BandPass_dict[key]['response'])
        x_filter_mask = (y_filter_array>np.max(y_filter_array)*0.75) #<TODO># threshold 0.75
        x_highlights = [np.min(x_filter_array[x_filter_mask]), np.max(x_filter_array[x_filter_mask])]
        #print('x_highlights =', x_highlights)
        ax.fill_between(x_highlights, [ax.get_ylim()[0],ax.get_ylim()[0]], [ax.get_ylim()[1],ax.get_ylim()[1]], color='gold', alpha=0.5)



def test_1():
    # 
    z = 3.0
    SED_dict = get_SED_Radio_mm_FIR_lines(L_FIR = 5e12, linewidth = 500.0, starburst = 0.0, z = z)
    
    fig = plt.figure(figsize=(7.0,4.0))
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(SED_dict['lambda_obs_um']/(1.+z), SED_dict['Snu_obs_mJy'])
    ax1.set_xlabel(r'Observing Wavelength [$\mu\mathrm{m}$]')
    ax1.set_ylabel(r'Flux Density [$\mathrm{mJy}$]')
    #ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim([30,3e5])
    ax1.set_ylim([0.01,1000])
    
    plot_Herschel_BandPass()
    
    #ax1.set_xlabel(r'$\lambda$ [$\mu m$]', fontsize=15)
    #ax1.set_ylabel(r'$S_{\nu}$ [$mJy$]', fontsize=15)
    #ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(x),0)))).format(x) if np.log10(x)<=3 else ''%(x))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.text(0.04, 0.90, r'$z=%0.2f$'%(z), transform=ax1.transAxes, fontsize=14)
    #ax1.text(0.04, 0.82, r'$\log \, M_\mathrm{d}=%0.2f$'%(np.log10(M_dust)), transform=ax1.transAxes, fontsize=14)
    #fig.tight_layout()
    fig.subplots_adjust(bottom=0.16, top=0.96, left=0.18, right=0.96)
    #plt.show(block=True)
    fig.savefig('Plot_test_FFT_1.pdf')
    os.system('open Plot_test_FFT_1.pdf')
    
    sys.exit()




if __name__ == '__main__':
    # 
    # 
    test_1()
    
    current_time_str = datetime.today().strftime('%Y-%m-%d %Hh%Mm%Ss %Z')
    
    list_linewidth = np.array([300.0, 500.0, 800.0])
    list_starburst = np.array([0.0, 1.0])
    
    meta = {}
    meta['NVAR'] = '2'
    meta['TVAR1'] = 'band wavelengths (rest-frame) (lambda) [um]'
    meta['TVAR2'] = 'S_nu_Jy * dL_Mpc**2 / M_dust. Final M_dust = a * dL**2 / (1.+z).'
    meta['CVAR1'] = '1 # the colomn number of VAR1'
    meta['CVAR2'] = '2 # the colomn number of VAR2'
    meta['NVAR1'] = 'NaN # NAXIS1'
    meta['NVAR2'] = '%d # NAXIS2'%(len(list_linewidth)*len(list_starburst))
    meta['NPAR'] = '3 # Parameter 1 2 3 ...'
    meta['TPAR1'] = 'linewidth'
    meta['TPAR2'] = 'starburst'
    meta['TPAR3'] = 'L_dust # non-independent. Final L_IR = TPAR3 * a * 4*pi*dL**2 / (1.+z).'
    meta['NPAR1'] = '%d # number of linewidth Values'%(len(list_linewidth))
    meta['NPAR2'] = '%d # number of starburst values'%(len(list_starburst))
    meta['NPAR3'] = '1 # non-independent'
    meta['CPAR1'] = '3 # column number of this parameter'
    meta['CPAR2'] = '4 # column number of this parameter'
    meta['CPAR3'] = '5 # column number of this parameter'
    
    lambda_um = []
    Snu_mJy = []
    linewidth = []
    starburst = []
    L_dust = []
    for i in range(len(list_linewidth)):
        for j in range(len(list_starburst)):
            print('list_linewidth[%d/%d], list_starburst[%d/%d]'%(i+1,len(list_linewidth),j+1,len(list_starburst)))
            SED = get_SED_Radio_mm_FIR_lines(linewidth = list_linewidth[i], starburst = list_starburst[j], verbose = False)
            lambda_um.extend(SED['lambda_um'])
            Snu_mJy.extend(SED['Snu_mJy'])
            linewidth.extend(SED['Snu_mJy']*0.0 + list_linewidth[i])
            starburst.extend(SED['Snu_mJy']*0.0 + list_starburst[j])
            inte_val = 0.0
            for k in range(1,len(SED['lambda_um'])):
                if SED['lambda_um'][k-1] >= 8.0 and SED['lambda_um'][k] <= 1000.0:
                    freq_step = np.abs(2.99792458e5/(SED['lambda_um'][k])-2.99792458e5/(SED['lambda_um'][k-1])) # GHz
                    inte_step = (SED['Snu_mJy'][k]+SED['Snu_mJy'][k-1])/2.0 * freq_step / 40.31970 # 1 Lsun Mpc-2 = 40.31970 mJy GHz
                    inte_val += inte_step
            L_dust.extend(SED['Snu_mJy']*0.0 + inte_val)
    # 
    meta['NVAR1'] = '%d # NAXIS1'%(len(SED['lambda_um']))
    # 
    otb = Table()
    otb['lambda_um'] = lambda_um
    otb['Snu_mJy'] = Snu_mJy
    otb['linewidth'] = linewidth
    otb['starburst'] = starburst
    otb['L_dust'] = L_dust
    
    otb['lambda_um'].format = '%0.6e'
    otb['Snu_mJy'].format = '%0.6e'
    otb['linewidth'].format = '%0.1f'
    otb['starburst'].format = '%0.1f'
    otb['L_dust'].format = '%0.6e'
    
    pprint(meta)
    #otb.meta = meta
    #with open('lib.MBB.SED.json', 'w') as fp:
    #    json.dump(meta, fp, indent=4, sort_keys=False)
    
    otb.write('lib.Radio.mm.FIR.lines.SED', format='ascii.fixed_width', delimiter=' ', overwrite=True)
    
    with open('lib.Radio.mm.FIR.lines.SED', 'r+') as fp:
        lines = fp.readlines() # read old content
        fp.seek(0) # go back to the beginning of the file
        fp.write('# %s\n'%(current_time_str))
        fp.write('# \n')
        for key in meta:
            fp.write('# %-6s = %s\n'%(key, meta[key])) # write new content at the beginning
        # 
        for i in range(len(lines)): # write old content after new
            if i == 0:
                fp.write('# \n')
                fp.write('# %s'%(re.sub(r'^  ', r'', lines[i])))
            else:
                fp.write(lines[i])

        
    
    print('Output to "lib.Radio.mm.FIR.lines.SED"!')
























