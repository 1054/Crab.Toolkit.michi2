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
from scipy.interpolate import interp1d
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = '13' # https://matplotlib.org/users/customizing.html
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



from makelibSED_Radio_mm_FIR_Lines import (get_SED_Radio_mm_FIR_lines, plot_Herschel_BandPass)


def test_1():
    # 
    z = 3.0
    SED_dict = get_SED_Radio_mm_FIR_lines(L_IR = 5e12, linewidth = 300.0, starburst = 0.0, z = z)
    
    grid = plt.GridSpec(4, 4) # NY = 4, NX = 4, 
    
    fig = plt.figure(figsize=(7.0,6.0))
    #ax1 = fig.add_subplot(2,1,1)
    ax1 = fig.add_subplot(grid[0:2,0:])
    ax1.plot(SED_dict['lambda_obs_um'], SED_dict['Snu_obs_mJy'])
    ax1.plot(SED_dict['lambda_obs_um'], SED_dict['Snu_obs_mJy']+0.0002)
    ax1.set_xlabel(r'Observing Wavelength [$\mu\mathrm{m}$]')
    ax1.set_ylabel(r'Flux Density [$\mathrm{mJy}$]')
    #ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim([30,3e5])
    ax1.set_ylim([0.0001,1000])
    plot_Herschel_BandPass(ax1)
    #ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(x),0)))).format(x) if np.log10(x)<=3 else ''%(x))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.text(0.04, 0.90, r'$z=%0.2f$'%(z), transform=ax1.transAxes, fontsize=14)
    
    freq_rest = 2.99792458e5/(SED_dict['lambda_obs_um'])*(1.+z)
    Snu_obs = SED_dict['Snu_obs_mJy']
    freq_rest = freq_rest[::-1]
    Snu_obs = Snu_obs[::-1]
    
    ax2 = fig.add_subplot(grid[2,0])
    ax2.plot(freq_rest, Snu_obs)
    ax2.plot(freq_rest, Snu_obs+0.0002)
    ax2.set_ylabel(r'Flux Density [$\mathrm{mJy}$]')
    ax2.set_xscale('log')
    ax2.set_xlim(np.array([115.2712018-2.0,115.2712018+2.2]))
    ax2.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax2.get_xlim()[0]), freq_rest < (ax2.get_xlim()[1]) ) ])])
    ax2.set_xticks(np.arange(115.2712018-2.0,115.2712018+2.2,1.0))
    ax2.text(0.04, 0.95, r'CO(1-0)', va='top', transform=ax2.transAxes, fontsize=12)
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax2.tick_params(axis='x', rotation=40)
    
    ax3 = fig.add_subplot(grid[2,1])
    ax3.plot(freq_rest, Snu_obs)
    ax3.plot(freq_rest, Snu_obs+0.0002)
    ax3.set_xscale('log')
    ax3.set_xlim(np.array([230.538-2.0,230.538+2.2]))
    ax3.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax3.get_xlim()[0]), freq_rest < (ax3.get_xlim()[1]) ) ])])
    ax3.set_xticks(np.arange(230.538-2.0,230.538+2.2,1.0))
    ax3.text(0.04, 0.95, r'CO(2-1)', va='top', transform=ax3.transAxes, fontsize=12)
    ax3.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax3.tick_params(axis='x', rotation=40)
    
    ax4 = fig.add_subplot(grid[2,2])
    ax4.plot(freq_rest, Snu_obs)
    ax4.plot(freq_rest, Snu_obs+0.0002)
    ax4.set_xscale('log')
    ax4.set_xlim(np.array([806.6518060-4.0,806.6518060+4.2]))
    ax4.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax4.get_xlim()[0]), freq_rest < (ax4.get_xlim()[1]) ) ])])
    ax4.set_xticks(np.arange(806.6518060-4.0,806.6518060+4.2,2.0))
    ax4.text(0.04, 0.95, r'CO(7-6)+[CI](2-1)', va='top', transform=ax4.transAxes, fontsize=11)
    ax4.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax4.tick_params(axis='x', rotation=40)
    
    ax5 = fig.add_subplot(grid[2,3])
    ax5.plot(freq_rest, Snu_obs)
    ax5.plot(freq_rest, Snu_obs+0.0002)
    ax5.set_xscale('log')
    ax5.set_xlim(np.array([1151.9854520-4.0,1151.9854520+4.2]))
    ax5.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax5.get_xlim()[0]), freq_rest < (ax5.get_xlim()[1]) ) ])])
    ax5.set_xticks(np.arange(1151.9854520-4.0,1151.9854520+4.2,4.0))
    ax5.text(0.04, 0.95, r'[CO](10-9)+H2O', va='top', transform=ax5.transAxes, fontsize=12)
    ax5.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax5.tick_params(axis='x', rotation=40)
    
    linefreq = 2459.38010; linelabel = r'[NII]122'
    linefreq = 1461.13141; linelabel = r'[NII]205'
    ax6 = fig.add_subplot(grid[3,0])
    ax6.plot(freq_rest, Snu_obs)
    ax6.plot(freq_rest, Snu_obs+0.0002)
    ax6.set_xscale('log')
    ax6.set_xlim(np.array([linefreq-8.0,linefreq+8.2]))
    ax6.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax6.get_xlim()[0]), freq_rest < (ax6.get_xlim()[1]) ) ])])
    ax6.set_xticks(np.arange(linefreq-8.0,linefreq+8.2,4.0))
    ax6.text(0.04, 0.95, linelabel, va='top', transform=ax6.transAxes, fontsize=12)
    ax6.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax6.tick_params(axis='x', rotation=40)
    
    linefreq = 1900.53690; linelabel = r'[CII]158'
    ax7 = fig.add_subplot(grid[3,1])
    ax7.plot(freq_rest, Snu_obs)
    ax7.plot(freq_rest, Snu_obs+0.0002)
    ax7.set_xscale('log')
    ax7.set_xlim(np.array([linefreq-8.0,linefreq+8.2]))
    ax7.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax7.get_xlim()[0]), freq_rest < (ax7.get_xlim()[1]) ) ])])
    ax7.set_xticks(np.arange(linefreq-8.0,linefreq+8.2,4.0))
    ax7.text(0.04, 0.95, linelabel, va='top', transform=ax7.transAxes, fontsize=12)
    ax7.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax7.tick_params(axis='x', rotation=40)
    
    linefreq = 3393.00624; linelabel = r'[OIII]88'
    #linefreq = 5785.87959; linelabel = r'[OIII]51' # or [OIII]52
    ax8 = fig.add_subplot(grid[3,2])
    ax8.plot(freq_rest, Snu_obs)
    ax8.plot(freq_rest, Snu_obs+0.0002)
    ax8.set_xscale('log')
    ax8.set_xlim(np.array([linefreq-8.0,linefreq+8.2]))
    ax8.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax8.get_xlim()[0]), freq_rest < (ax8.get_xlim()[1]) ) ])])
    ax8.set_xticks(np.arange(linefreq-8.0,linefreq+8.2,4.0))
    ax8.text(0.04, 0.95, linelabel, va='top', transform=ax8.transAxes, fontsize=12)
    ax8.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax8.tick_params(axis='x', rotation=40)
    
    linefreq = 2.99792458e5/63.183705; linelabel = r'[OI]63'
    ax9 = fig.add_subplot(grid[3,3])
    ax9.plot(freq_rest, Snu_obs)
    ax9.plot(freq_rest, Snu_obs+0.0002)
    ax9.set_xscale('log')
    ax9.set_xlim(np.array([linefreq-8.0,linefreq+8.2]))
    ax9.set_ylim([0.0, 1.8 * np.max(Snu_obs[ np.logical_and(freq_rest > (ax9.get_xlim()[0]), freq_rest < (ax9.get_xlim()[1]) ) ])])
    ax9.set_xticks(np.arange(linefreq-8.0,linefreq+8.2,4.0))
    ax9.text(0.04, 0.95, linelabel, va='top', transform=ax9.transAxes, fontsize=12)
    ax9.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('%0.3f'%(x)) )  ) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax9.tick_params(axis='x', rotation=40)
    
    plt.tight_layout(rect=[0.0,0.045,1.0,1.0]) # (left, bottom, right, top) 
    
    fig.text(0.5, 0.035, r'Rest-frame Frequency [$\mathrm{GHz}$]', ha='center', fontsize=13)
    
    #fig.subplots_adjust(bottom=0.16, top=0.96, left=0.18, right=0.96)
    #plt.show(block=True)
    fig.savefig('Plot_test_1.pdf')
    os.system('open Plot_test_1.pdf')
    
    sys.exit()




if __name__ == '__main__':
    # 
    # 
    test_1()
























