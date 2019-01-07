#!/usr/bin/env python
# 
# generate radio/mm/FIR strong lines as template SEDs, incl. CO, [CI], [CII], [OIII], [OI], [NII]
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
from scipy.interpolate import interp1d
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
mpl.rcParams['text.usetex'] = True

sys.path.append(os.path.expanduser('~')+os.sep+'Cloud/Github/Crab.Toolkit.Python/lib/crab/crabpdbi/')
from CrabPdBI import ( find_radio_lines_in_frequency_range, 
                       calc_radio_line_flux_from_IR_luminosity, 
                       convert_flux2lprm, 
                       convert_flux2lsun, 
                     )



global_wavelength_um_log10_resolution = 0.01
global_velocity_resolution = 10.0 #<TODO> 



def get_SED_Radio_mm_FIR_lines(L_IR, linewidth = 500.0, starburst = 0.0, z = 0.0, verbose = True):
    if verbose:
        print('get_SED_Radio_mm_FIR_lines(L_IR = %e, linewidth = %.1f, starburst = %.1f, z = %.3f)'%(L_IR, linewidth, starburst, z))
    
    if z <= 0.0:
        raise ValueError('Error! z is invalid!')
        
    # 
    # get all lines
    output_line_freqs, output_line_names = find_radio_lines_in_frequency_range([2.99792458e5/1e5, 2.99792458e5/30.0], Redshift = 0.0, set_output_line_names = True, include_faint_lines = False)
    if verbose:
        print('get_SED_Radio_mm_FIR_lines() found following lines:')
        for i in np.argsort(output_line_freqs)[::-1]:
            print('    %-20s at %18.8f GHz'%(output_line_names[i], output_line_freqs[i]))
    output_line_specs = []
    output_wavelength_um = []
    output_Snu_mJy = []
    for i in np.argsort(output_line_freqs)[::-1]:
        # 
        IR_color = 0.6 + starburst #<TODO># IR_color
        IR_luminosity = L_IR
        line_name = output_line_names[i]
        rest_freq = output_line_freqs[i]
        obs_freq = rest_freq / (1.+z)
        # 
        velocity_resolution = global_velocity_resolution
        velocity_grid = np.arange(+2500.0, -2500.0-velocity_resolution, -velocity_resolution) #<TODO> 
        frequency_grid = velocity_grid/2.99792458e5 * obs_freq + obs_freq
        wavelength_grid = 2.99792458e5/frequency_grid # in increasing direction
        # 
        dL = cosmo.luminosity_distance(z).to(u.Mpc).value
        line_flux = calc_radio_line_flux_from_IR_luminosity(line_name, IR_luminosity, z, starburstiness = starburst, IR_color = IR_color, verbose = False)
        line_lprm = convert_flux2lprm(line_flux, rest_freq, z)
        # 
        line_sigma = linewidth / (2.0*np.sqrt(2.0*np.log(2.0))) # km/s, from Gaussian FWHM to Gaussian sigma
        line_spec = 1.0/(line_sigma*np.sqrt(2.0*np.pi)) * np.exp(-(velocity_grid)**2/(2.0*line_sigma**2))
        line_spec = line_spec * line_flux*1e3 # mJy
        # 
        output_line_specs.append({'obs_wavelength_um':wavelength_grid.tolist(), 
                                  'obs_Sv_mJy':line_spec.tolist(), 
                                  'obs_velocity_grid':velocity_grid.tolist(), 
                                  'obs_frequency_grid':frequency_grid.tolist(), 
                                  'interpolator':interp1d(wavelength_grid, line_spec), 
                                  'line_name':line_name,
                                  'rest_freq':rest_freq,
                                })
        # print debug info
        if verbose:
            print('Line name %s, wavelength range %.3f - %.3f um, line lprm %e K km/s pc2, L_IR %e L_sun, line flux %s Jy km/s.'%(line_name, wavelength_grid[0], wavelength_grid[-1], line_lprm, IR_luminosity, line_flux))
        # merge wavelength grid, remove intersection
        if len(output_wavelength_um) == 0:
            output_wavelength_um = wavelength_grid.tolist()
            output_Snu_mJy = line_spec.tolist()
        else:
            if verbose:
                print('Line name %s, output wavelength array range %.3f - %.3f um.'%(line_name, output_wavelength_um[0], output_wavelength_um[-1]))
            # extend with all wavelength_grid elements which are not close to already existing elements in output_wavelength_um
            # see -- https://stackoverflow.com/questions/32513424/find-intersection-of-numpy-float-arrays
            intersecting_atol = np.abs(2.99792458e5/obs_freq - 2.99792458e5/((velocity_resolution)/2.99792458e5*obs_freq+obs_freq)) / 2.0
            non_intersected_index = ~np.isclose(    np.array(output_wavelength_um)[:,None],\
                                                    wavelength_grid,\
                                                    atol=intersecting_atol\
                                               ).any(0)
            #print(list(zip(wavelength_grid,non_intersected_index))) #<DEBUG>
            output_wavelength_um.extend(wavelength_grid[non_intersected_index].tolist())
            output_Snu_mJy.extend(line_spec[non_intersected_index].tolist())
            if verbose:
                print('Line name %s, len(output_wavelength_um) = %d, non_intersected_index = %d, atol = %s'%(line_name, len(output_wavelength_um), np.count_nonzero(non_intersected_index), intersecting_atol ))
            # 
            # if there has intersection, we need to coadd the spectrum
            if np.count_nonzero(~non_intersected_index) > 0:
                x = np.array(output_wavelength_um)
                y = np.array(output_Snu_mJy)
                x2 = wavelength_grid[~non_intersected_index]
                y2 = line_spec[~non_intersected_index]
                mask = np.logical_and(x>=np.min(x2), x<=np.max(x2))
                y[mask] += interp1d(x2,y2)(x[mask])
                output_wavelength_um = x.tolist()
                output_Snu_mJy = y.tolist()
                if verbose:
                    print('Line name %s, coadding intersected wavelength range %.3f - %.3f um.'%(line_name, x2[0], x2[-1]))
            
    # 
    # add broad uniform wavelength grid points
    broad_wavelength_um_log10_resolution = global_wavelength_um_log10_resolution
    broad_wavelength_um_log10_grid = np.arange(np.log10(30.0),np.log10(1e5)+broad_wavelength_um_log10_resolution,broad_wavelength_um_log10_resolution)
    intersecting_atol = broad_wavelength_um_log10_resolution / 100.0 # np.log10(np.abs(2.99792458e5/obs_freq - 2.99792458e5/((global_velocity_resolution/2.0)/2.99792458e5*obs_freq+obs_freq)))
    non_intersected_index = ~np.isclose(    np.log10(np.array(output_wavelength_um))[:,None],\
                                            broad_wavelength_um_log10_grid,\
                                            atol=intersecting_atol\
                                       ).any(0) # index for broad_wavelength_um_log10_grid
    output_wavelength_um.extend(np.power(10.0,broad_wavelength_um_log10_grid[non_intersected_index]).tolist())
    output_Snu_mJy.extend([0.0]*len(broad_wavelength_um_log10_grid[non_intersected_index])) # combine broad and fine wavelength grid
    output_wavelength_um, output_Snu_mJy = zip(*sorted(zip(output_wavelength_um, output_Snu_mJy))) # sort by wavelength
    
    wavelength_um_log10 = np.log10(np.array(output_wavelength_um))
    wavelength_um = np.array(output_wavelength_um)
    nu = 2.99792458e5/wavelength_um
    Snu = np.array(output_Snu_mJy)
    
    # 
    # loop again to spline line fluxes
    #for i in range(len(output_line_specs)):
    #    line_mask = np.logical_and(loop_wavelength_um >= np.min(output_line_specs[i]['obs_wavelength_um']), 
    #                               loop_wavelength_um < np.max(output_line_specs[i]['obs_wavelength_um']) )
    #    Snu[line_mask] += output_line_specs[i]['interpolator'](wavelength_um[line_mask])
    
    # 
    # 
    SED_dict = {}
    SED_dict['lambda_um'] = wavelength_um / (1.+z)
    SED_dict['Snu_mJy'] = Snu / (1.+z)
    SED_dict['Snu_Jy'] = Snu/1e3 / (1.+z)
    SED_dict['lambda_rest_um'] = wavelength_um / (1.+z)
    SED_dict['Snu_rest_mJy'] = Snu / (1.+z)
    SED_dict['Snu_rest_Jy'] = Snu/1e3 / (1.+z)
    if z > 0.0:
        SED_dict['lambda_obs_um'] = wavelength_um
        SED_dict['Snu_obs_mJy'] = Snu
        SED_dict['Snu_obs_Jy'] = Snu/1e3
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
                BandPass_dict[key]['lambda_um'].append(float(t[0])) # already converted AA to um in these files
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
    SED_dict = get_SED_Radio_mm_FIR_lines(L_IR = 5e12, linewidth = 500.0, starburst = 0.0, z = z)
    
    fig = plt.figure(figsize=(7.0,4.0))
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(SED_dict['lambda_obs_um'], SED_dict['Snu_obs_mJy'])
    ax1.set_xlabel(r'Observing Wavelength [$\mu\mathrm{m}$]')
    ax1.set_ylabel(r'Flux Density [$\mathrm{mJy}$]')
    #ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim([30,3e5])
    ax1.set_ylim([0.0001,1000])
    
    plot_Herschel_BandPass(ax1)
    
    #ax1.set_xlabel(r'$\lambda$ [$\mu m$]', fontsize=15)
    #ax1.set_ylabel(r'$S_{\nu}$ [$mJy$]', fontsize=15)
    #ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(x),0)))).format(x) if np.log10(x)<=3 else ''%(x))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.text(0.04, 0.90, r'$z=%0.2f$'%(z), transform=ax1.transAxes, fontsize=14)
    #ax1.text(0.04, 0.82, r'$\log \, M_\mathrm{d}=%0.2f$'%(np.log10(M_dust)), transform=ax1.transAxes, fontsize=14)
    #fig.tight_layout()
    fig.subplots_adjust(bottom=0.16, top=0.96, left=0.18, right=0.96)
    #plt.show(block=True)
    fig.savefig('Plot_test_1.pdf')
    os.system('open Plot_test_1.pdf')
    
    sys.exit()




if __name__ == '__main__':
    # 
    # 
    #test_1()
    
    current_time_str = datetime.today().strftime('%Y-%m-%d %Hh%Mm%Ss %Z')
    
    list_linewidth = np.array([300.0, 500.0, 800.0])
    list_starburst = np.array([0.0, 1.0])
    
    meta = {}
    meta['NVAR'] = '2'
    meta['TVAR1'] = 'band wavelengths (rest-frame) (lambda) [um]'
    meta['TVAR2'] = 'S_nu_mJy for L_IR=1e12 at z=3 and removed (1.+z) dependency. Note that L^{\\Prime}_line = 3.25e7 \\int{S_nu_obs d nu} * dL^2 * nu_rest^{-2} / (1.+z) (Solomon 2005), so fitted model * PAR4^2 / dL^2 * (1.+z) to the get physical line spectrum.'
    meta['CVAR1'] = '1 # the colomn number of VAR1'
    meta['CVAR2'] = '2 # the colomn number of VAR2'
    meta['NVAR1'] = 'NaN # NAXIS1'
    meta['NVAR2'] = '%d # NAXIS2'%(len(list_linewidth)*len(list_starburst))
    meta['NPAR'] = '5 # Parameter 1 2 3 ...'
    meta['TPAR1'] = 'linewidth'
    meta['TPAR2'] = 'starburst'
    meta['TPAR3'] = 'z # redshift used for generating the model, non-independent.'
    meta['TPAR4'] = 'dL # luminosity distance used for generating the model, non-independent.'
    meta['TPAR5'] = 'LIR # IR luminosity used for generating the model, non-independent. Final L_IR = a * PAR4^2 / dL^2 * (1.+z) * PAR5, where a is the fitted normalization of this model, and z and dL are the redshift and luminosity distance of the input source.'
    meta['NPAR1'] = '%d # number of linewidth Values'%(len(list_linewidth))
    meta['NPAR2'] = '%d # number of starburst values'%(len(list_starburst))
    meta['NPAR3'] = '1 # non-independent'
    meta['NPAR4'] = '1 # non-independent'
    meta['NPAR5'] = '1 # non-independent'
    meta['CPAR1'] = '3 # column number of this parameter'
    meta['CPAR2'] = '4 # column number of this parameter'
    meta['CPAR3'] = '5 # column number of this parameter'
    meta['CPAR4'] = '6 # column number of this parameter'
    meta['CPAR5'] = '7 # column number of this parameter'
    
    lambda_um = []
    Snu_mJy = []
    linewidth = []
    starburst = []
    listLIR = []
    listz = []
    listdL = []
    for i in range(len(list_linewidth)):
        for j in range(len(list_starburst)):
            print('list_linewidth[%d/%d], list_starburst[%d/%d]'%(i+1,len(list_linewidth),j+1,len(list_starburst)))
            
            L_IR = 1e12
            z = 3.0
            dL = cosmo.luminosity_distance(z).to(u.Mpc).value
            SED = get_SED_Radio_mm_FIR_lines(L_IR = L_IR, linewidth = list_linewidth[i], starburst = list_starburst[j], z = z, verbose = False)
            
            lambda_um.extend(SED['lambda_um'])
            Snu_mJy.extend(SED['Snu_mJy'])
            linewidth.extend(SED['Snu_mJy']*0.0 + list_linewidth[i])
            starburst.extend(SED['Snu_mJy']*0.0 + list_starburst[j])
            listz.extend(SED['Snu_mJy']*0.0 + z)
            listdL.extend(SED['Snu_mJy']*0.0 + dL)
            listLIR.extend(SED['Snu_mJy']*0.0 + (L_IR))
    # 
    meta['NVAR1'] = '%d # NAXIS1'%(len(SED['lambda_um']))
    # 
    otb = Table()
    otb['lambda_um'] = lambda_um
    otb['Snu_mJy'] = Snu_mJy
    otb['linewidth'] = linewidth
    otb['starburst'] = starburst
    otb['z'] = listz
    otb['dL'] = listdL
    otb['LIR'] = listLIR
    
    otb['lambda_um'].format = '%0.6e'
    otb['Snu_mJy'].format = '%0.6e'
    otb['linewidth'].format = '%0.1f'
    otb['starburst'].format = '%0.1f'
    otb['z'].format = '%0.3f'
    otb['dL'].format = '%0.3f'
    otb['LIR'].format = '%0.6e'
    
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
























