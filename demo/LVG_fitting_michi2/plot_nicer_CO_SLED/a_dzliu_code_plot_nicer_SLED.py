#!/usr/bin/env python
# 
from __future__ import print_function
import os, sys, re, json, copy, datetime
from astropy.table import Table
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
#from scipy.interpolate import spline # spline was removed from version 1.3.1
from scipy.interpolate import InterpolatedUnivariateSpline

def spline(x_model, y_model, x_input, order = 3):
    #print('x_model, y_model', x_model, y_model)
    spmodel = InterpolatedUnivariateSpline(x_model, y_model, k = order+1) # k=3 mean cubic
    return spmodel(x_input)

def map_X_species_to_CO_J(X_species):
    if np.isscalar(X_species):
        return map_X_species_to_CO_J([X_species])[0]
    else:
        CO_J = []
        for i in range(len(X_species)):
            if X_species[i] > 101000000 and X_species[i] < 102000000:
                CO_J.append(int(np.round((X_species[i]-101000000.)/1000.)))
            else:
                raise Exception('Error! Unrecognized X_species %d for the CO line!'%(X_species[i]))
        CO_J = np.array(CO_J)
        return CO_J

def map_X_species_to_CI_J(X_species):
    if np.isscalar(X_species):
        return map_X_species_to_CI_J([X_species])[0]
    else:
        CI_J = []
        for i in range(len(X_species)):
            if np.isclose(X_species[i], 102001000, rtol = 1e-10, atol = 0.1):
                CI_J.append(1)
            elif np.isclose(X_species[i], 102002001, rtol = 1e-10, atol = 0.1):
                CI_J.append(2)
            elif np.isclose(X_species[i], 102002000, rtol = 1e-10, atol = 0.1):
                CI_J.append(3)
            else:
                raise Exception('Error! Unrecognized X_species %d for the [CI] line!'%(X_species[i]))
        CI_J = np.array(CI_J)
        return CI_J


# 
# set common renorm CO J
common_renorm_CO_J = 2
# 
# prepare figure
fig = plt.figure(figsize=(5.2,3.8))
plt.subplots_adjust(left=0.11, right=0.98, bottom=0.21, top=0.97)
ax = fig.add_subplot(1,1,1)
# 
# read tables
tb0 = Table.read('flux_co_ci.txt', format='ascii')
tb1 = Table.read('result_two_component_fit/dump/1/LVG_LIB1', format='ascii')
tb2 = Table.read('result_two_component_fit/dump/1/LVG_LIB2', format='ascii')
OBS_CO_mask = np.logical_and(tb0.columns[0] > 101000000., tb0.columns[0] < 102000000.)
OBS_CI_mask = np.logical_and(tb0.columns[0] > 102000000., tb0.columns[0] < 103000000.)
LIB1_CO_mask = np.logical_and(tb1.columns[0] > 101000000., tb1.columns[0] < 102000000.)
LIB1_CI_mask = np.logical_and(tb1.columns[0] > 102000000., tb1.columns[0] < 103000000.)
LIB2_CO_mask = np.logical_and(tb2.columns[0] > 101000000., tb2.columns[0] < 102000000.)
LIB2_CI_mask = np.logical_and(tb2.columns[0] > 102000000., tb2.columns[0] < 103000000.)
LIB1_CO_J = map_X_species_to_CO_J(tb1.columns[0][LIB1_CO_mask]) # ((tb1.columns[0][LIB1_CO_mask]-101000000.)/1000.).astype(int)
LIB1_CI_J = map_X_species_to_CI_J(tb1.columns[0][LIB1_CI_mask]) # np.arange(np.count_nonzero(LIB1_CI_mask)).astype(int) + 1 # ((tb1.columns[0][LIB1_CI_mask]-102000000.)/1000.).astype(int)[0:2]
LIB1_CO_F = tb1.columns[1][LIB1_CO_mask]
LIB1_CI_F = tb1.columns[1][LIB1_CI_mask]
LIB2_CO_J = map_X_species_to_CO_J(tb2.columns[0][LIB2_CO_mask]) # ((tb2.columns[0][LIB2_CO_mask]-101000000.)/1000.).astype(int)
LIB2_CI_J = map_X_species_to_CI_J(tb2.columns[0][LIB2_CI_mask]) # np.arange(np.count_nonzero(LIB2_CI_mask)).astype(int) + 1 # ((tb2.columns[0][LIB2_CI_mask]-102000000.)/1000.).astype(int)[0:2]
LIB2_CO_F = tb2.columns[1][LIB2_CO_mask]
LIB2_CI_F = tb2.columns[1][LIB2_CI_mask]
#<DEBUG>#for i in range(len(LIB1_CO_J)):
#<DEBUG>#    print(LIB1_CO_J[i], LIB1_CO_F[i], LIB2_CO_F[i])
# 
# read observed data points
OBS_CO_J = map_X_species_to_CO_J(tb0.columns[0][OBS_CO_mask])
OBS_CO_F = tb0.columns[1][OBS_CO_mask].data
OBS_CO_FErr = tb0.columns[2][OBS_CO_mask].data
OBS_CI_J = map_X_species_to_CI_J(tb0.columns[0][OBS_CI_mask])
OBS_CI_F = tb0.columns[1][OBS_CI_mask].data
OBS_CI_FErr = tb0.columns[2][OBS_CI_mask].data
# 
# renormalize observed data points
print('OBS_CO_J', OBS_CO_J)
print('OBS_CO_F', OBS_CO_F)
OBS_renorm = 1.0 / OBS_CO_F[OBS_CO_J==common_renorm_CO_J]
print('OBS_renorm', OBS_renorm)
OBS_CO_F = OBS_CO_F*OBS_renorm
OBS_CO_FErr = OBS_CO_FErr*OBS_renorm
OBS_CI_F = OBS_CI_F*OBS_renorm
OBS_CI_FErr = OBS_CI_FErr*OBS_renorm
# 
# spline model grid
LIB_CO_J = np.linspace(1.0, 10.0, num=50, endpoint=True)
LIB_CO_F = spline(LIB1_CO_J, LIB1_CO_F+LIB2_CO_F, LIB_CO_J)
OFF_CI_J = np.max(LIB_CO_J) # axis offset
LIB_CI_J = np.array([1, 2])
LIB_CI_F = spline(LIB1_CI_J, LIB1_CI_F+LIB2_CI_F, LIB_CI_J, order=0)
# 
# renormalize model curves
LIB_renorm = OBS_renorm # 1.0 / spline(LIB_CO_J, LIB_CO_F, common_renorm_CO_J)
LIB_CO_F = LIB_CO_F * LIB_renorm
LIB_CI_F = LIB_CI_F * LIB_renorm
LIB1_CO_F = LIB1_CO_F * LIB_renorm
LIB2_CO_F = LIB2_CO_F * LIB_renorm
LIB1_CI_F = LIB1_CI_F * LIB_renorm
LIB2_CI_F = LIB2_CI_F * LIB_renorm
# 
# plot model curves
ax.plot(LIB_CO_J, spline(LIB1_CO_J, LIB1_CO_F, LIB_CO_J), color='blue', lw=0.9, label='Best-fit low-exc.')
ax.plot(LIB_CO_J, spline(LIB2_CO_J, LIB2_CO_F, LIB_CO_J), color='red', lw=0.9, label='Best-fit high-exc.')
ax.plot(LIB_CO_J, LIB_CO_F, color='k', lw=1.8, label='Best-fit total')
ax.plot(OFF_CI_J+LIB_CI_J, spline(LIB1_CI_J, LIB1_CI_F, LIB_CI_J, order=0), color='blue', lw=0.9, label='__none__')
ax.plot(OFF_CI_J+LIB_CI_J, spline(LIB2_CI_J, LIB2_CI_F, LIB_CI_J, order=0), color='red', lw=0.9, label='__none__')
ax.plot(OFF_CI_J+LIB_CI_J, LIB_CI_F, color='k', lw=1.8, label='__none__')
#<DEBUG>#
#ax.plot(LIB1_CO_J, LIB1_CO_F+LIB2_CO_F, color='k', lw=1.8, marker='x', label='__none__')
# 
# plot observed data points
#ax.errorbar(OBS_CO_J, OBS_CO_F, OBS_CO_FErr, linestyle='', marker='o', markersize=10, markerfacecolor='none', markeredgecolor='k', color='k', markeredgewidth=1.8, lw=1.8, elinewidth=1.8, capsize=5, 
#            label='ID 2299 Narrow \n(This work)', zorder=30)
mask = (OBS_CO_F>3.0*OBS_CO_FErr)
if np.count_nonzero(mask) > 0:
    ax.errorbar(OBS_CO_J[mask], OBS_CO_F[mask], OBS_CO_FErr[mask], linestyle='', marker='o', markersize=10, markerfacecolor='none', markeredgecolor='k', color='k', markeredgewidth=1.8, lw=1.8, elinewidth=1.8, capsize=5, 
                label='Target\n(This work)', zorder=30)
if np.count_nonzero(~mask) > 0:
    ax.errorbar(OBS_CO_J[~mask], 3.0*OBS_CO_FErr[~mask], OBS_CO_FErr[~mask], uplims=True, color='k', elinewidth=1.8, capthick=1.8, capsize=5, 
                label='__none__', zorder=30)
mask = (OBS_CI_F>3.0*OBS_CI_FErr)
#ax.errorbar(OFF_CI_J+OBS_CI_J, OBS_CI_F, OBS_CI_FErr, linestyle='', marker='o', markersize=10, markerfacecolor='none', markeredgecolor='k', color='k', markeredgewidth=1.8, lw=1.8, elinewidth=1.8, capsize=5, 
#            label='__none__', zorder=30)
if np.count_nonzero(mask) > 0:
    ax.errorbar(OFF_CI_J+OBS_CI_J[mask], OBS_CI_F[mask], OBS_CI_FErr[mask], linestyle='', marker='o', markersize=10, markerfacecolor='none', markeredgecolor='k', color='k', markeredgewidth=1.8, lw=1.8, elinewidth=1.8, capsize=5, 
                label='__none__', zorder=30)
if np.count_nonzero(~mask) > 0:
    ax.errorbar(OFF_CI_J+OBS_CI_J[~mask], 3.0*OBS_CI_FErr[~mask], OBS_CI_FErr[~mask], uplims=True, color='k', elinewidth=1.8, capthick=1.8, capsize=5, 
                label='__none__', zorder=30)

# 
# plot local SB CO SLED
#tbX = Table.read('datatable_CO_SLED_local_SB_Daddi2015.txt', format='ascii')
#tbX_renorm = spline(LIB_CO_J, LIB_CO_F, tbX['avLine'][0]) / tbX['avFlux'][0]
#ax.fill_between(tbX['avLine'], (tbX['avFlux']-tbX['avFErr'])*tbX_renorm, (tbX['avFlux']+tbX['avFErr'])*tbX_renorm, color='orange', alpha=0.5)
# 
# plot HFLS-3
tbX = Table.read('datatable_CO_SLED_SMG_HFLS3_Riechers2013.txt', format='ascii')
#tbX_renorm = spline(LIB_CO_J, LIB_CO_F, tbX['j'][0]) / tbX['f'][0]
tbX_renorm = 1.0 / tbX['f'][tbX['j']==common_renorm_CO_J]
ax.plot(tbX['j'], tbX['f']*tbX_renorm, dashes=(3,1), marker='x', mew=1.8, color='orange', alpha=0.6, label=r'HFLS-3')
# 
# plot HFLS-3 ([CI])
tbX = Table.read('datatable_CI_SLED_SMG_HFLS3_Riechers2013.txt', format='ascii')
#tbX_renorm = spline(LIB_CI_J, LIB_CI_F, tbX['j'][0], order=0) / tbX['f'][0] # using CO renorm
ax.plot(OFF_CI_J+tbX['j'], tbX['f']*tbX_renorm, dashes=(3,1), marker='x', mew=1.8, color='gold', alpha=0.6, label=r'__none__')
# 
# plot BzK
tbX = Table.read('datatable_CO_SLED_BzK_Daddi2015.txt', format='ascii')
#tbX_renorm = spline(LIB_CO_J, LIB_CO_F, tbX['j'][0]) / tbX['f'][0]
tbX_renorm = 1.0 / tbX['f'][tbX['j']==common_renorm_CO_J]
ax.plot(tbX['j'], tbX['f']*tbX_renorm, dashes=(3,1), linestyle='none', color='#509e00', marker='*', markersize=12, mew=1.8, mec='none', mfc='#509e00', alpha=0.6, label=r'BzKs')
# 
# plot Milky Way
tbX = Table.read('datatable_CO_SLED_MilkyWay_Fixsen1999.txt', format='ascii')
#tbX_renorm = spline(LIB_CO_J, LIB_CO_F, tbX['j'][0]) / tbX['f'][0]
tbX_renorm = 1.0 / tbX['f'][tbX['j']==common_renorm_CO_J]
ax.plot(tbX['j'], tbX['f']*tbX_renorm, dashes=(2,2), linestyle='none', color='#0000cc', marker='^', markersize=7, mew=1.8, mec='none', mfc='#0000cc', alpha=0.6, label=r'MW')
# 
# set axis interval
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
ax.grid(True, c='#cccccc', ls='dotted', lw=0.5)
# 
# show legend
ax.legend(loc='upper right', fontsize=10.0, framealpha=0.7, borderaxespad=0.2, borderpad=0.2, handletextpad=0.15, labelspacing=0.3, frameon=True)
# 
# xtitle
ax.set_xlim([0.05, np.max(OFF_CI_J+LIB_CI_J)+0.5])
ax.set_ylabel(r'Relative Line Flux', fontsize=13, labelpad=10) # $\mathrm{[Jy\,km\,s^{-1}]}$
ax.tick_params(which='both', labelsize=12)
ax.tick_params(axis='x', rotation=60, labelsize=11.5)
def axis_tick_func_formatter(x, pos):
    if x <= 10:
        return 'CO(%d-%d)'%(x,x-1)
    else:
        return '[CI](%d-%d)'%(x-10,x-10-1)
ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axis_tick_func_formatter))
# 
# ylog
#ax.set_yscale('log')
ax.set_ylim([-0.19, 5.5])
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1.0))
# 
# tight_layout
#fig.tight_layout()
# 
#plt.show()
fig.savefig('Plot_nicer_SLED.pdf', transparent=True, dpi=300)
os.system('open Plot_nicer_SLED.pdf')





