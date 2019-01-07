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
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from copy import copy
from pprint import pprint
#from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)
from astropy import units as u
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

from makelibSED_MBB import get_SED_MBB, get_SED_CMB

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Radio_mm_FIR_Lines')

from makelibSED_Radio_mm_FIR_Lines import get_Herschel_BandPass, plot_Herschel_BandPass, get_SED_Radio_mm_FIR_lines



def get_Filter_Curves():
    Filter_Curves_dict = {}
    key = 'AzTEC'
    Filter_Curves_dict[key] = {}
    Filter_Curves_dict[key]['freq_GHz'] = np.linspace(240.0, 290.0, num=11, endpoint=True) # see Figure 4 of https://arxiv.org/pdf/0801.2783.pdf
    Filter_Curves_dict[key]['lambda_um'] = 2.99792458e5/Filter_Curves_dict[key]['freq_GHz']
    Filter_Curves_dict[key]['response'] = Filter_Curves_dict[key]['freq_GHz']*0.0 + 1.0 #<TODO># assuming flat response on flux density in units of Jy
    key = '1.2mm'
    Filter_Curves_dict[key] = {}
    Filter_Curves_dict[key]['freq_GHz'] = np.linspace(2.99792458e5/1200.0-5.0, 2.99792458e5/1200.0+5.0, num=11, endpoint=True) # assuming 10 GHz
    Filter_Curves_dict[key]['lambda_um'] = 2.99792458e5/Filter_Curves_dict[key]['freq_GHz']
    Filter_Curves_dict[key]['response'] = Filter_Curves_dict[key]['freq_GHz']*0.0 + 1.0 #<TODO># assuming flat response on flux density in units of Jy
    key = '1.3mm'
    Filter_Curves_dict[key] = {}
    Filter_Curves_dict[key]['freq_GHz'] = np.linspace(2.99792458e5/1300.0-5.0, 2.99792458e5/1300.0+5.0, num=11, endpoint=True) # assuming 10 GHz
    Filter_Curves_dict[key]['lambda_um'] = 2.99792458e5/Filter_Curves_dict[key]['freq_GHz']
    Filter_Curves_dict[key]['response'] = Filter_Curves_dict[key]['freq_GHz']*0.0 + 1.0 #<TODO># assuming flat response on flux density in units of Jy
    key = '3mm'
    Filter_Curves_dict[key] = {}
    Filter_Curves_dict[key]['freq_GHz'] = np.linspace(2.99792458e5/3000.0-5.0, 2.99792458e5/3000.0+5.0, num=11, endpoint=True) # assuming 10 GHz
    Filter_Curves_dict[key]['lambda_um'] = 2.99792458e5/Filter_Curves_dict[key]['freq_GHz']
    Filter_Curves_dict[key]['response'] = Filter_Curves_dict[key]['freq_GHz']*0.0 + 1.0 #<TODO># assuming flat response on flux density in units of Jy
    key = 'SCUBA2_850'
    Filter_Curves_dict[key] = {}
    Filter_Curves_dict[key]['freq_GHz'] = np.linspace(2.99792458e5/850.0-5.0, 2.99792458e5/850.0+5.0, num=11, endpoint=True) # assuming 10 GHz
    Filter_Curves_dict[key]['lambda_um'] = 2.99792458e5/Filter_Curves_dict[key]['freq_GHz']
    Filter_Curves_dict[key]['response'] = Filter_Curves_dict[key]['freq_GHz']*0.0 + 1.0 #<TODO># assuming flat response on flux density in units of Jy
    # 
    Herschel_BandPass_dict = get_Herschel_BandPass()
    for key in Herschel_BandPass_dict:
        Filter_Curves_dict[key] = Herschel_BandPass_dict[key]
    # 
    return Filter_Curves_dict


def test():
    z = 5.85
    M_dust = 10**8.0
    beta = 2.2
    T_dust = 40.0 # this is the observed T_dust, not phyiscal, nor intrinsic.
    T_dust_physical = 41.0 # take what Georgios has fitted including CMB
    L_IR = 8e11
    dL = cosmo.luminosity_distance(z).value
    
    tb_real_data = Table.read('/Users/dzliu/Work/DeepFields/Works_cooperated/2018_Shuowen_Jin/20181001_photometry_from_Shuowen/ID1929_z5.85.txt', format='ascii.commented_header')
    
    # compute CMB SED
    SED_CMB = get_SED_CMB(z, beta = beta, M_dust = M_dust)
    x_CMB = SED_CMB['lambda_obs_um']
    y_CMB = SED_CMB['Snu_obs_mJy']
    
    # compute physical dust SED if were no CMB
    SED_dust_physical = get_SED_MBB(beta = beta, T_dust = T_dust_physical, z = z, M_dust = M_dust, verbose = True)
    x_dust_physical = SED_dust_physical['lambda_obs_um']
    y_dust_physical = SED_dust_physical['Snu_obs_mJy']
    
    # compute observed flux, i.e., physical - CMB
    x_dust_against_CMB = y_dust_physical
    y_dust_against_CMB = y_dust_physical - y_CMB
    
    # compute intrinsic dust temperature if no CMB
    T_dust_intrinsic = np.power(T_dust_physical**(4+beta) - (cosmo.Tcmb0.value*(1.+z))**(4+beta), 1.0/(4+beta))
    SED_dust_intrinsic = get_SED_MBB(beta = beta, T_dust = T_dust_intrinsic, z = z, M_dust = M_dust, verbose = True)
    x_dust_intrinsic = SED_dust_intrinsic['lambda_obs_um']
    y_dust_intrinsic = SED_dust_intrinsic['Snu_obs_mJy']
    
    # compute diff
    x_diff_intrinsic_to_apparent = x_dust_intrinsic
    y_diff_intrinsic_to_apparent = y_dust_intrinsic - y_dust_against_CMB
    
    #intersecting_rtol = (x_dust_physical[1]-x_dust_physical[0])/20.0
    #non_intersected_index = ~np.isclose(    x_dust_physical[:,None],\
    #                                        x_CO,\
    #                                        rtol=intersecting_rtol\
    #                                   ).any(0)
    #x_dust_physical_plus_CO = x_dust_physical.tolist()
    #x_dust_physical_plus_CO.extend(x_CO[non_intersected_index].tolist())
    #y_dust_physical_plus_CO = y_dust_physical.tolist()
    #y_dust_physical_plus_CO.extend(y_CO[non_intersected_index].tolist()) # original CO line width add dust SED
    #x_dust_physical_plus_CO, y_dust_physical_plus_CO = zip(*sorted(zip(x_dust_physical_plus_CO, y_dust_physical_plus_CO))) # sort by wavelength
    
    SED_CO = get_SED_Radio_mm_FIR_lines(L_IR = L_IR, linewidth = 500.0, starburst = 0.0, z = z)
    x_CO = SED_CO['lambda_obs_um'] # already redshifted in get_SED_Radio_mm_FIR_lines()
    y_CO = SED_CO['Snu_obs_mJy'] # already redshifted in get_SED_Radio_mm_FIR_lines()
    x_dust_plus_CO_physical = x_CO
    y_dust_plus_CO_physical = y_CO + interp1d(x_dust_physical,y_dust_physical)(x_dust_plus_CO_physical)
    
    x_dust_plus_CO_against_CMB = x_dust_plus_CO_physical
    y_dust_plus_CO_against_CMB = y_dust_plus_CO_physical - interp1d(x_CMB,y_CMB)(x_dust_plus_CO_against_CMB)
    
    fig = plt.figure(figsize=(7.0,4.0))
    ax1 = fig.add_subplot(1,1,1)
    
    ax1.plot(x_CMB, 
             y_CMB, ls='dashed', lw=1.0, c='blue', alpha=0.8, label=r'CMB at $z=%0.3f$'%(z))
    
    ax1.plot(x_dust_intrinsic, 
             y_dust_intrinsic, ls='dashed', lw=2.0, c='darkgray', alpha=0.8, label=r'dust MBB ($\beta=%0.1f$, $T_{\mathrm{d}}^{\mathrm{in.}}=%0.1f$)'%(beta, T_dust_intrinsic))
    
    ax1.plot(x_dust_physical, 
             y_dust_physical, ls='dotted', lw=2.0, c='darkgray', alpha=0.8, label=r'dust MBB ($\beta=%0.1f$, $T_{\mathrm{d}}=%0.1f$)'%(beta, T_dust_physical))
    
    ax1.plot(x_dust_plus_CO_against_CMB, 
             y_dust_plus_CO_against_CMB, ls='solid', lw=1.0, c='k', alpha=0.8, label=r'dust - CMB + CO (obs.)')
    
    ax1.plot(x_diff_intrinsic_to_apparent, 
             y_diff_intrinsic_to_apparent, ls='dashed', lw=1.0, c='#00FF00', alpha=0.8, label=r'diff. intrins.-obs.')
    
    data_mask = (tb_real_data['flux'] > 3.0 * tb_real_data['dflux'])
    ax1.errorbar(tb_real_data['wave'][data_mask], tb_real_data['flux'][data_mask], tb_real_data['dflux'][data_mask], ls='none', c='magenta', alpha=0.7, marker='s', ms=8, capsize=6, capthick=1.2, label='obs. data')
    ax1.errorbar(tb_real_data['wave'][~data_mask], 3.0*tb_real_data['dflux'][~data_mask], yerr=0.5*tb_real_data['dflux'][~data_mask], uplims=True, ls='none', mfc='none', c='magenta', ms=8, capsize=6, capthick=1.2, label='__none__')
    
    # for each real data point, compute CO through filter curve as a continuum
    Filter_Curves_dict = get_Filter_Curves()
    for i in range(len(tb_real_data)):
        for key in Filter_Curves_dict:
            filter_lambda_um = np.array(Filter_Curves_dict[key]['lambda_um'])
            filter_response = np.array(Filter_Curves_dict[key]['response'])
            filter_range_mask = (filter_response >= np.max(filter_response)*0.75)
            obs_data_wavelength = tb_real_data['wave'][i]
            obs_data_filter_mask = np.logical_and(obs_data_wavelength>=np.min(filter_lambda_um[filter_range_mask]), obs_data_wavelength<=np.max(filter_lambda_um[filter_range_mask]))
            if np.count_nonzero(obs_data_filter_mask) > 0:
                # check whether x_CO covers the most part of filter curve
                x_CO_does_cover_filter = np.logical_and(np.min(x_CO)<=np.min(filter_lambda_um[filter_range_mask]), np.max(x_CO)>=np.max(filter_lambda_um[filter_range_mask]))
                if x_CO_does_cover_filter:
                    print('computing flux through filter %s for obs data at wavelength %s um'%(key, obs_data_wavelength))
                    x_CO_filter_mask = np.logical_and(x_CO>=np.min(filter_lambda_um), x_CO<=np.max(filter_lambda_um))
                    x_dust_filter_mask = np.logical_and(x_dust_against_CMB>=np.min(filter_lambda_um), x_dust_against_CMB<=np.max(filter_lambda_um))
                    if np.count_nonzero(x_CO_filter_mask) > 0:
                        # 
                        # compute for low resolution dust+CMB SED
                        filter_response_interp = (interp1d(filter_lambda_um, filter_response))(x_dust_against_CMB[x_dust_filter_mask])
                        y_dust_against_CMB_through_filter = np.sum(y_dust_against_CMB[x_dust_filter_mask] * filter_response_interp) / \
                                                            np.sum(filter_response_interp)
                        # 
                        # compute for high resolution CO+dust+CMB SED
                        filter_response_interp = (interp1d(filter_lambda_um, filter_response))(x_CO[x_CO_filter_mask])
                        y_dust_plus_CO_against_CMB_through_filter = np.sum(y_dust_plus_CO_against_CMB[x_CO_filter_mask] * filter_response_interp) / \
                                                                    np.sum(filter_response_interp)
                        # 
                        # compute CO contribution to the continuum
                        y_CO_as_a_continuum = y_dust_plus_CO_against_CMB_through_filter - y_dust_against_CMB_through_filter
                        # 
                        # print debug info
                        print('computing flux through filter %s for obs data at wavelength %s um: %s mJy (CO as a continuum contributes %s mJy)'%(key, obs_data_wavelength, y_dust_plus_CO_against_CMB_through_filter, y_CO_as_a_continuum))
                        # 
                        # plot flux through filter
                        ax1.plot([np.min(filter_lambda_um[filter_range_mask]),np.max(filter_lambda_um[filter_range_mask])], 
                                 [y_dust_plus_CO_against_CMB_through_filter,y_dust_plus_CO_against_CMB_through_filter], 
                                 ls='solid', lw=5.0, c='k', alpha=0.5, label=r'__none__', zorder=10)
                        #break
    
    #plot_Herschel_BandPass(ax1)
    
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim([1e2,2e4])
    #ax1.set_xlim([1000,1300])
    ax1.set_ylim([1e-2,1e2])
    ax1.set_xlabel(r'$\lambda$ [$\mu m$]', fontsize=15)
    ax1.set_ylabel(r'$S_{\nu}$ [$mJy$]', fontsize=15)
    ax1.set_yticks(np.power(10.0,np.arange(-2,3,1)))
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(x),0)))).format(x) if np.abs(np.log10(x))<=3 else '$10^{%d}$'%(np.log10(x)) ) ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y) if np.abs(np.log10(y))<=3 else '$10^{%d}$'%(np.log10(y)) ) ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.tick_params(labelsize=12)
    ax1.tick_params(direction='in', axis='both', which='both')
    ax1.tick_params(top=True, right=True, which='both') # top='on', right='on' -- deprecated -- use True/False instead
    ax1.grid(True, ls='--', lw=0.25)
    ax1.text(0.04, 0.90, r'$z=%0.2f$'%(z), transform=ax1.transAxes, fontsize=14)
    ax1.text(0.04, 0.82, r'$\log \, M_\mathrm{d}=%0.2f$'%(np.log10(M_dust)), transform=ax1.transAxes, fontsize=14)
    fig.subplots_adjust(bottom=0.16, top=0.96, left=0.18, right=0.96)
    #plt.show(block=True)
    fig.savefig('Plot_MBB_test_4_CMB_CO.pdf')
    os.system('open Plot_MBB_test_4_CMB_CO.pdf')
    
    sys.exit()




if __name__ == '__main__':
    # 
    # 
    test()
























