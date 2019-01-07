#!/usr/bin/env python
# 
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
    z = 10.0
    M_dust = 10**8.0
    beta = 2.0
    T_dust_intrinsic = 18.0
    T_dust = np.power(T_dust_intrinsic**(4+beta) + (cosmo.Tcmb0.value*(1.+z))**(4+beta), 1.0/(4+beta))
    
    SED_dust = get_SED_MBB(beta = beta, T_dust = T_dust, M_dust = M_dust, z = z, verbose = True)
    SED_dust_intrinsic = get_SED_MBB(beta = beta, T_dust = T_dust_intrinsic, M_dust = M_dust, z = z, verbose = True)
    SED_CMB = get_SED_CMB(z)
    # 
    fig = plt.figure(figsize=(7.0,4.0))
    ax1 = fig.add_subplot(1,1,1)
    dL = cosmo.luminosity_distance(z).value
    x_dust_physical = SED_dust['lambda_um'] * (1.+z)
    y_dust_physical = SED_dust['Snu_mJy'] * M_dust / dL**2 * (1.+z)
    x_dust_intrinsic = SED_dust_intrinsic['lambda_um'] * (1.+z)
    y_dust_intrinsic = SED_dust_intrinsic['Snu_mJy'] * M_dust / dL**2 * (1.+z)
    x_CMB = SED_CMB['lambda_um'] # already redshifted in get_SED_CMB()
    y_CMB = SED_CMB['Snu_mJy'] * M_dust / dL**2 # already redshifted in get_SED_CMB()
    x_dust_against_CMB = x_dust_physical
    y_dust_against_CMB = y_dust_physical - y_CMB
    
    x_diff_intrinsic_to_apparent = x_dust_against_CMB
    y_diff_intrinsic_to_apparent = y_dust_intrinsic - y_dust_against_CMB
    
    
    ax1.plot(x_CMB/(1.+z), 
             y_CMB, ls='dashed', lw=2.0, c='blue', alpha=0.8, label=r'CMB at $z=%0.3f$'%(z))
    
    ax1.plot(x_dust_intrinsic/(1.+z), 
             y_dust_intrinsic, ls='solid', lw=2.0, c='darkgray', alpha=0.8, label=r'dust MBB ($\beta=%0.1f$, $T_{\mathrm{d}}^{\mathrm{in}}=%0.1f$)'%(beta, T_dust_intrinsic))
    
    ax1.plot(x_dust_physical/(1.+z), 
             y_dust_physical, ls='solid', lw=1.0, c='k', alpha=0.8, label=r'dust MBB ($\beta=%0.1f$, $T_{\mathrm{d}}=%0.1f$)'%(beta, T_dust))
    
    ax1.plot(x_dust_against_CMB/(1.+z), 
             y_dust_against_CMB, ls='dashed', lw=2.0, c='darkgray', alpha=0.8, label=r'dust against CMB, i.e., obs.')
    
    ax1.plot(x_diff_intrinsic_to_apparent/(1.+z), 
             y_diff_intrinsic_to_apparent, ls='dashed', lw=1.0, c='#00FF00', alpha=0.8, label=r'diff. in.-appa.')
    
    #plot_Herschel_BandPass(ax1)
    
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim([1e1,1e6])
    #ax1.set_xlim([1000,1300])
    ax1.set_ylim([1e-8,1e2])
    ax1.set_xlabel(r'$\lambda_{\mathrm{rest}}$ [$\mu m$]', fontsize=15)
    ax1.set_ylabel(r'$S_{\nu}$ [$mJy$]', fontsize=15)
    ax1.set_yticks(np.power(10.0,np.arange(-2,2,1)))
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
    fig.savefig('Plot_MBB_test_5_CMB_z10.pdf')
    os.system('open Plot_MBB_test_5_CMB_z10.pdf')
    
    sys.exit()




if __name__ == '__main__':
    # 
    # 
    test()
























