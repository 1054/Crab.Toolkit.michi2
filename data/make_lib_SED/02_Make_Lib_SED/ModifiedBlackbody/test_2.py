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

from makelibSED_MBB import get_SED_MBB




def test():
    z = 5.85
    M_dust = 10**8.0
    beta = 2.2
    T_dust = 41.0
    SED = get_SED_MBB(beta = beta, T_dust = T_dust, M_dust = M_dust, z = z, verbose = True)
    SED_ref = get_SED_MBB(beta = 1.8, T_dust = 25.0, M_dust = M_dust, z = z, verbose = True)
    SED_comp1 = get_SED_MBB(beta = 2.8, T_dust = 25.0, M_dust = M_dust, z = z, verbose = True)
    SED_comp2 = get_SED_MBB(beta = 1.8, T_dust = 50.0, M_dust = M_dust, z = z, verbose = True)
    # 
    tb_real_data = Table.read('/Users/dzliu/Work/DeepFields/Works_cooperated/2018_Shuowen_Jin/20181001_photometry_from_Shuowen/ID1929_z5.85.txt', format='ascii.commented_header')
    # 
    fig = plt.figure(figsize=(7.0,4.0))
    ax1 = fig.add_subplot(1,1,1)
    dL = cosmo.luminosity_distance(z).value
    x = SED['lambda_um'] * (1.+z)
    y = SED['Snu_Jy'] * M_dust / dL**2 * (1.+z) * 1e3
    x_ref = SED_ref['lambda_um'] * (1.+z)
    y_ref = SED_ref['Snu_Jy'] * M_dust / dL**2 * (1.+z) * 1e3
    x_comp1 = SED_comp1['lambda_um'] * (1.+z)
    y_comp1 = SED_comp1['Snu_Jy'] * M_dust / dL**2 * (1.+z) * 1e3
    x_comp2 = SED_comp2['lambda_um'] * (1.+z)
    y_comp2 = SED_comp2['Snu_Jy'] * M_dust / dL**2 * (1.+z) * 1e3
    ax1.plot(x_ref, y_ref, ls='solid', lw=1.0, c='k', alpha=0.8, label=r'$\beta=1.8$, $T_{\mathrm{d}}=25.0$')
    ax1.plot(x_comp1, y_comp1, ls='dashed', lw=1.0, c='blue', alpha=0.8, label=r'$\beta=2.8$, $T_{\mathrm{d}}=25.0$')
    ax1.plot(x_comp2, y_comp2, ls='dashed', lw=1.0, c='red', alpha=0.8, label=r'$\beta=1.8$, $T_{\mathrm{d}}=50.0$')
    ax1.plot(x, y, ls='dotted', lw=2.0, c='darkgray', alpha=0.8, label=r'$\beta=%0.1f$, $T_{\mathrm{d}}=%0.1f$'%(beta, T_dust))
    data_mask = (tb_real_data['flux'] > 3.0 * tb_real_data['dflux'])
    #ax1.plot(tb_real_data['wave'][data_mask], tb_real_data['flux'][data_mask], ls='none', mfc='magenta', mec='none', marker='s', ms=8, alpha=0.5, label='__none__')
    ax1.errorbar(tb_real_data['wave'][data_mask], tb_real_data['flux'][data_mask], tb_real_data['dflux'][data_mask], ls='none', c='magenta', alpha=0.7, marker='s', ms=8, capsize=6, capthick=1.2, label='obs. data')
    ax1.errorbar(tb_real_data['wave'][~data_mask], 3.0*tb_real_data['dflux'][~data_mask], yerr=0.5*tb_real_data['dflux'][~data_mask], uplims=True, ls='none', mfc='none', c='magenta', ms=8, capsize=6, capthick=1.2, label='__none__')
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim([1e2,1e4])
    ax1.set_ylim([1e-2,1e2])
    ax1.set_xlabel(r'$\lambda$ [$\mu m$]', fontsize=15)
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
    fig.savefig('Plot_MBB_test_2.pdf')
    os.system('open Plot_MBB_test_2.pdf')
    
    sys.exit()




if __name__ == '__main__':
    # 
    # 
    test()
























