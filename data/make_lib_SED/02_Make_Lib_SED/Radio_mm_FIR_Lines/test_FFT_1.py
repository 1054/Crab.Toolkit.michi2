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
from CrabPdBI import ( calc_radio_line_flux_from_IR_luminosity, 
                       )








def test_FFT_1():
    # 
    #x = 1 + 2 * 1j
    dXWaveNumber = 0.1 * u.cm**-1
    XWaveNumber = np.arange(1.0,100.0,dXWaveNumber.value) * u.cm**-1 # cm-1
    XWavelength = (1.0 / XWaveNumber).to(u.um) # um
    XFrequency = 2.99792458e5/XWavelength.value * u.GHz # GHz
    dXFrequency = 2.99792458e5/dXWaveNumber.value * u.GHz # GHz
    #print(XWaveNumber[0], ',', XWavelength[0])
    
    array1 = np.zeros(XWaveNumber.shape,dtype=np.complex_)
    array1_FT = np.zeros(XWaveNumber.shape,dtype=np.complex_)
    
    # insert some lines in data domain
    idx = np.abs(XFrequency-115.0*u.GHz).argmin()
    if idx >= 0:
        array1[idx] = 1.0
    #
    #array1_FT = fft(array1)
    for i in range(len(XWaveNumber)):
        array1_FT[i] = np.sum(array1 * np.exp(-2*np.pi*1j*(XWaveNumber[i].to(u.cm**-1).value)*(XFrequency.to(u.GHz).value)) * dXFrequency.to(u.GHz).value)
    
    # add some feature to FT domain
    array1_FT = array1_FT * 0.0
    #array1_FT.real[ np.logical_and(XWaveNumber>10/u.cm, XWaveNumber<50/u.cm) ] = 1.0
    array1_FT.real = array1_FT.real + 10.0
    #array1_FT = (array1_FT + 1.0) * np.exp(-2*np.pi*1j * ((XWaveNumber.value*10.0)))
    #array1_FT = (array1_FT + 1.0) * np.cos(2*np.pi*theta) + 1j * np.sin(2*np.pi*theta)
    # 
    #array1 = ifft(array1_FT)
    for i in range(len(XFrequency)):
        array1[i] = np.sum(array1_FT * np.exp(+2*np.pi*1j*(XWaveNumber.to(u.cm**-1).value)*(XFrequency[i].to(u.GHz).value)) * dXWaveNumber.to(u.cm**-1).value)
    
    fig = plt.figure(figsize=(7.0,8.0))
    ax1 = fig.add_subplot(4,2,1)
    
    ax1.plot(XWaveNumber, array1_FT.real)
    ax1.text(0.04, 0.82, r'FT real', transform=ax1.transAxes, fontsize=14)
    
    ax3 = fig.add_subplot(4,2,3)
    ax3.plot(XWaveNumber, array1_FT.imag)
    ax3.text(0.04, 0.82, r'FT imag', transform=ax3.transAxes, fontsize=14)
    
    ax5 = fig.add_subplot(4,2,5)
    ax5.plot(XWaveNumber, np.abs(array1_FT))
    ax5.text(0.04, 0.82, r'FT ampl', transform=ax5.transAxes, fontsize=14)
    
    ax7 = fig.add_subplot(4,2,7)
    ax7.plot(XWaveNumber, np.angle(array1_FT))
    ax7.text(0.04, 0.82, r'FT phase', transform=ax7.transAxes, fontsize=14)
    ax7.set_xlabel(r'wavenumber [$\mathrm{cm^{-1}}$]')
    
    ax2 = fig.add_subplot(4,2,2)
    ax2.plot(XFrequency, array1.real)
    ax2.text(0.04, 0.82, r'data real', transform=ax2.transAxes, fontsize=14)
    
    ax4 = fig.add_subplot(4,2,4)
    ax4.plot(XFrequency, array1.imag)
    ax4.text(0.04, 0.82, r'data imag', transform=ax4.transAxes, fontsize=14)
    
    ax6 = fig.add_subplot(4,2,6)
    ax6.plot(XFrequency, np.abs(array1))
    ax6.text(0.04, 0.82, r'data ampl', transform=ax6.transAxes, fontsize=14)
    
    ax8 = fig.add_subplot(4,2,8)
    ax8.plot(XFrequency, np.angle(array1))
    ax8.text(0.04, 0.82, r'data phase', transform=ax8.transAxes, fontsize=14)
    ax8.set_xlabel(r'frequency [$\mathrm{GHz}$]')
    
    #ax1.legend()
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    #ax1.set_xlim([0.1,3e5])
    #ax1.set_ylim([1e-5,1e4])
    #ax1.set_xlabel(r'$\lambda$ [$\mu m$]', fontsize=15)
    #ax1.set_ylabel(r'$S_{\nu}$ [$mJy$]', fontsize=15)
    #ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(x),0)))).format(x) if np.log10(x)<=3 else ''%(x))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    #ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    #ax1.tick_params(labelsize=12)
    #ax1.tick_params(direction='in', axis='both', which='both')
    #ax1.tick_params(top=True, right=True, which='both') # top='on', right='on' -- deprecated -- use True/False instead
    #ax1.grid(True, ls='--', lw=0.25)
    #ax1.text(0.04, 0.90, r'$z=%0.2f$'%(z), transform=ax1.transAxes, fontsize=14)
    #ax1.text(0.04, 0.82, r'$\log \, M_\mathrm{d}=%0.2f$'%(np.log10(M_dust)), transform=ax1.transAxes, fontsize=14)
    fig.tight_layout()
    #fig.subplots_adjust(bottom=0.16, top=0.96, left=0.18, right=0.96)
    #plt.show(block=True)
    fig.savefig('Plot_test_FFT_1.pdf')
    os.system('open Plot_test_FFT_1.pdf')
    
    sys.exit()




if __name__ == '__main__':
    # 
    # 
    test_FFT_1()
    
























