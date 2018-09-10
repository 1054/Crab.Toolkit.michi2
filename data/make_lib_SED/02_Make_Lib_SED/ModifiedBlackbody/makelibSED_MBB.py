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



def get_SED_MBB(beta = 1.8, T_dust = 25.0, M_dust = 0.0, z = 0.0, verbose = False):
    # MBB: modified blackbody, modified by a power law in frequency, ν^β (Hildebrand 1983)
    # 
    # setup wavelength grid
    rest_wavelength_um_log10 = np.arange(np.log10(3.0), np.log10(1e6)+0.05, 0.05)
    rest_wavelength_um = np.power(10.0, rest_wavelength_um_log10)
    nu = 2.99792458e5/rest_wavelength_um # GHz
    # 
    # planck function blackbody
    rest_Bv = blackbody_nu(rest_wavelength_um * u.um, T_dust * u.K) # erg s^-1 cm^-2 Hz^-1 Sr^-1
    rest_Bv_at_100um = blackbody_nu(100.0 * u.um, T_dust * u.K) # erg s^-1 cm^-2 Hz^-1 Sr^-1, ~= 1e-12
    rest_Bv_at_850um = blackbody_nu(850.0 * u.um, T_dust * u.K) # erg s^-1 cm^-2 Hz^-1 Sr^-1
    rest_Bv_at_230GHz = blackbody_nu((2.99792458e5/230.0) * u.um, T_dust * u.K) # erg s^-1 cm^-2 Hz^-1 Sr^-1
    # 
    # opacity, kappa0 = 1-np.exp(-tau0), tau0 = -np.log(1-kappa0)
    #nu0 = 3000.0 # GHz # ν0 is the frequency where optical depth equals unity (Draine 2006), see https://arxiv.org/pdf/1206.1595.pdf Eq(1) and text below Eq(1). The theoretically expected value of ν0 is 3 THz (i.e. λ0=100 µm), although this value is unconstrained by data (see discussion in Conley et al. 2011). 
    nu0 = 2.99792458e5 / 850.0 # the reference frequency, 850um
    kappa0 = 1.5 * (u.cm**2/u.g)
                    # κ850 = 0.15 m^2 kg^-1 = 0.15 * 1e4 / 1e3 cm^2 g-1 (Weingartner & Draine 2001; Dunne et al. 2003). # κν is the dust mass absorption coefficient at ν. See https://arxiv.org/pdf/1206.1595.pdf Eq(8) and the text below Eq(8). 
                    # see also U et al. (2012) Eq(5)
                    # We can also see a temperature-depentend kappa -- http://www2.iap.fr/users/srinivas/Dust2Galaxies/Gordon.pdf -- Page 17
    #nu0 = 230.0 * u.GHz # GHz # the reference frequency, 230GHz, 1.3mm
    #kappa0 = 0.009 * (u.cm**2/u.g)
                    # κ230 = 0.009 cm^2/g (Ossenkopf & Henning 1994),  which accounts for the dust-to-gas ratio, so that the free parameter in the fit is the gas column density N. See https://arxiv.org/pdf/1203.0025.pdf Eq(1) and text below Eq(1). 
                    # κ230 = 0.009 cm^2 g^-1 is the emissivity at 230 GHz of the dust grains at a gas density of 10^6 cm−3 covered by a thin ice mantle (Ossenkopf & Henning 1994, Column 6 of Table 1). -- see https://arxiv.org/pdf/1408.5429.pdf Sect. 3. 
                    # 
    #n0 = 2.99792458e5 / 1000. * u.GHz
    #kappa0 = 0.77 * (u.cm**2/u.g)
                    # Values of kappa_nu at a conventional frequency of around 1 mm are in the range 0.04-0.15 m2 kg-1 (Hughes, 1996). 
                    # Dunne et al. (2000) adopt a value of 0.077 m2 kg-1. See -- https://ned.ipac.caltech.edu/level5/Sept04/Blain/Blain2_2.html, Sect. 2.2.1. 
                    # 
    # 
    # 
    # modified blackbody
    # tau == µH2 mH κν N(H2), where µ~2.8 (Kauffmann et al. 2008) is the mean weight of Mmol per mH. -- see https://arxiv.org/pdf/1101.4654.pdf
    # κν == κ0(ν/ν0)^β
    tau = kappa0.value * np.power(nu/nu0, beta)
    rest_Sv = rest_Bv * (1.0 - np.exp(-tau)) # erg s^-1 g^-1 Hz^-1 Sr^-1
    # 
    tau_at_100um = kappa0.value * np.power((2.99792458e5/100.0)/nu0, beta)
    rest_Sv_at_100um = rest_Bv_at_100um * (1.0 - np.exp(-tau_at_100um)) # erg s^-1 g^-1 Hz^-1 Sr^-1, ~= 1e-14
    # 
    tau_at_850um = kappa0.value * np.power((2.99792458e5/850.0)/nu0, beta)
    rest_Sv_at_850um = rest_Bv_at_850um * (1.0 - np.exp(-tau_at_850um)) # erg s^-1 g^-1 Hz^-1 Sr^-1
    # 
    tau_at_230GHz = kappa0.value * np.power((230.0)/nu0, beta)
    rest_Sv_at_230GHz = rest_Bv_at_230GHz * (1.0 - np.exp(-tau_at_230GHz)) # erg s^-1 g^-1 Hz^-1 Sr^-1
    # 
    # S_v = kappa_v * B_v(T) * M_dust / dL^2
    #dL = 1*u.Mpc # Mpc
    #dL2 = (dL.to(u.cm))**2# (dL)**2 * 9.52140e48 # Mpc^2 -> cm^2
    if verbose:
        print('rest_Bv_at_100um', rest_Bv_at_100um, '# erg s^-1 Hz^-1 sr^-1')
        print('rest_Bv_at_850um', rest_Bv_at_850um, '# erg s^-1 Hz^-1 sr^-1')
        print('rest_Bv_at_230GHz', rest_Bv_at_230GHz, '# erg s^-1 Hz^-1 sr^-1')
        print('tau_at_100um', tau_at_100um, '1-exp(-tau)', (1.0 - np.exp(-tau_at_100um)))
        print('tau_at_850um', tau_at_850um, '1-exp(-tau)', (1.0 - np.exp(-tau_at_850um)))
        print('tau_at_230GHz', tau_at_230GHz, '1-exp(-tau)', (1.0 - np.exp(-tau_at_230GHz)))
        print('rest_Sv_at_100um', rest_Sv_at_100um, '# erg s^-1 g^-1 Hz^-1')
        print('rest_Sv_at_850um', rest_Sv_at_850um, '# erg s^-1 g^-1 Hz^-1')
        print('rest_Sv_at_230GHz', rest_Sv_at_230GHz, '# erg s^-1 g^-1 Hz^-1')
        print('nu0', nu0, '# GHz')
        print('kappa0', kappa0, '# cm^2 g^-1')
    #print('dL', dL, '# Mpc')
    #print('dL2', dL2, '# cm^2')
    rest_Sv = rest_Sv / ((u.Mpc).to(u.cm)/u.Mpc)**2 * u.sr * ((u.solMass).to(u.g)/u.solMass) * 1e23*u.Jy/(u.erg/u.s/u.cm/u.cm/u.Hz)
    rest_Sv_at_100um = rest_Sv_at_100um / ((u.Mpc).to(u.cm)/u.Mpc)**2 * u.sr * ((u.solMass).to(u.g)/u.solMass) * 1e23*u.Jy/(u.erg/u.s/u.cm/u.cm/u.Hz)
    rest_Sv_at_850um = rest_Sv_at_850um / ((u.Mpc).to(u.cm)/u.Mpc)**2 * u.sr * ((u.solMass).to(u.g)/u.solMass) * 1e23*u.Jy/(u.erg/u.s/u.cm/u.cm/u.Hz)
    rest_Sv_at_230GHz = rest_Sv_at_230GHz / ((u.Mpc).to(u.cm)/u.Mpc)**2 * u.sr * ((u.solMass).to(u.g)/u.solMass) * 1e23*u.Jy/(u.erg/u.s/u.cm/u.cm/u.Hz)
    if verbose:
        print('rest_Sv_at_100um', rest_Sv_at_100um, '# Jy Msun^-1 Mpc^2')
        print('rest_Sv_at_850um', rest_Sv_at_850um, '# Jy Msun^-1 Mpc^2')
        print('rest_Sv_at_230GHz', rest_Sv_at_230GHz, '# Jy Msun^-1 Mpc^2')
    # 
    if M_dust > 0 and z > 0:
        # Msun
        dL = cosmo.luminosity_distance(z) # Mpc
        if verbose:
            print('dL', dL, '# Mpc')
        obs_Sv = rest_Sv * M_dust*u.solMass / dL**2
        obs_Sv_at_100um = rest_Sv_at_100um * M_dust*u.solMass / dL**2
        obs_Sv_at_850um = rest_Sv_at_850um * M_dust*u.solMass / dL**2
        obs_Sv_at_230GHz = rest_Sv_at_230GHz * M_dust*u.solMass / dL**2
        if verbose:
            print('obs_Sv_at_100um', obs_Sv_at_100um, '# Jy')
            print('obs_Sv_at_850um', obs_Sv_at_850um, '# Jy')
            print('obs_Sv_at_230GHz', obs_Sv_at_230GHz, '# Jy')
        
        obs_Sv = obs_Sv * (1.+z) # nu_obs * obs_Sv = nu_rest * rest_Sv, nu_rest = nu_obs * (1.+z)
        obs_Sv_at_100um = obs_Sv_at_100um * (1.+z)
        obs_Sv_at_850um = obs_Sv_at_850um * (1.+z)
        obs_Sv_at_230GHz = obs_Sv_at_230GHz * (1.+z)
        if verbose:
            print('obs_Sv_at_100um', obs_Sv_at_100um, '# Jy')
            print('obs_Sv_at_850um', obs_Sv_at_850um, '# Jy')
            print('obs_Sv_at_230GHz', obs_Sv_at_230GHz, '# Jy')
        
    # 
    # Or according to https://arxiv.org/pdf/1203.0025.pdf
    # S_v = Omega * N * kappa_0 * tau * B_v(T)
    # Omega = Area/dL^2 is the solid angle of the source
    # N is gas mass column density along the light of sight, thus Area * N = M_gas, thus S_v = Mdust/dL^2 * kappa_0 * tau * B_v(T)
    # κ230 = 0.009 cm^2 g^-1 is the emissivity at 230 GHz of the dust grains at a gas density of 106 cm−3 covered by a thin ice mantle (Ossenkopf & Henning 1994, Column 6 of Table 1). -- see https://arxiv.org/pdf/1408.5429.pdf Sect. 3. 
    # 
    # M_dust ~ 1e8 Msun ~ 2e41 g
    # Area ~ pi * 10**2 kpc^2 ~ 9e45 cm^2
    # dL^2 ~ 2e4**2 * 9e48 ~ 3.6e57 cm^2
    # Omega ~ 9e45/3.6e57 ~ 2.5e-12 sr
    # and
    # N ~ 2e21 cm^-2 ~ 4e-24 * 2e21 ~ 9e-3 g cm^-2
    # M_dust ~ N * Area / GDR ~ 2e41 g, consistent
    # and
    # kappa0 ~ 9e-3 cm^2 g^-1
    # B_v_230GHz ~ 8.182499051e-7 * 230**2 * 25 * 4.25e10 ~ 4.6e10 # Jy (230 GHz) (25 K) (1 sr = 4.25e10 arcsec^2)
    # so 
    # S_v ~ 2.5e-12 * 9e-3 * 9e-3 * 4.6e10 ~ 9e-6 Jy
    # so 
    # M_dust ~= N * Area / GDR = S_v / Omega / (kappa_0 * tau * B_v(T)) * Area / GDR = S_v / (kappa_0 * tau * B_v(T)) * dL^2 / GDR
    # 
    # In another method (Casey et al 2011; U et al. 2012; Casey et al. 2012;), M_dust ~= S_v * dL^2 / (kappa * B_v)
    # B_v ~ 3e-13 # [erg s-1 cm-2 Hz-1 sr-1] 
    # B_v = 2 * 6.62606957e-27 * 230e9**3 / 3e10**2 * (1.0 / ( np.exp( (6.62606957e-27*230e9)/(1.3806488e-16*25) ) - 1.0 )) # ~ 3e-13
    # S_v ~ 1 [mJy] ~ 1e-26 [erg s-1 cm-2 Hz-1]
    # dL^2 ~ (2e4)**2 * 9.52140e48 ~ 3.8e57 # [cm^2]
    # M_dust ~ S_v * dL^2 / (kappa * B_v) / 2e33 ~ 4e10 # [Msun] --> WRONG?
    # 
    SED_dict = {}
    SED_dict['lambda_um'] = rest_wavelength_um
    SED_dict['Snu_Jy'] = rest_Sv.value
    return SED_dict


def test():
    z = 5.05
    M_dust = 1e9
    beta = 1.8
    T_dust = 25.0
    SED = get_SED_MBB(beta = beta, T_dust = T_dust, M_dust = M_dust, z = z, verbose = True)
    # 
    fig = plt.figure(figsize=(7.0,4.0))
    ax1 = fig.add_subplot(1,1,1)
    dL = cosmo.luminosity_distance(z).value
    x = SED['lambda_um'] * (1.+z)
    y = SED['Snu_Jy'] * M_dust / dL**2 * (1.+z) * 1e3
    ax1.plot(x, y)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim([1e-6,1e3])
    ax1.set_xlabel(r'$\lambda$ [$\mu m$]', fontsize=15)
    ax1.set_ylabel(r'$S_{\nu}$ [$mJy$]', fontsize=15)
    #ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(x),0)))).format(x) if np.log10(x)<=3 else ''%(x))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    #ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y))) # https://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
    ax1.tick_params(labelsize=12)
    ax1.tick_params(direction='in', axis='both', which='both')
    ax1.tick_params(top=True, right=True, which='both') # top='on', right='on' -- deprecated -- use True/False instead
    ax1.grid(True, ls='--', lw=0.25)
    ax1.text(0.04, 0.90, '$z=%0.2f$'%(z), transform=ax1.transAxes, fontsize=14)
    ax1.text(0.04, 0.82, '$M_{dust}=%0.1e$'%(M_dust), transform=ax1.transAxes, fontsize=14)
    fig.subplots_adjust(bottom=0.16, top=0.96, left=0.18, right=0.96)
    plt.show(block=True)




if __name__ == '__main__':
    # 
    # 
    #test()
    
    list_beta = np.arange(1.5, 2.5+0.1, 0.1)
    list_T_dust = np.arange(10.0, 100.0+2.0, 2.0)
    
    meta = {}
    meta['NVAR'] = '2'
    meta['TVAR1'] = 'band wavelengths (rest-frame) (lambda) [um]'
    meta['TVAR2'] = 'S_nu_Jy * dL_Mpc**2 / M_dust. Final M_dust = a * dL**2 / (1.+z).'
    meta['CVAR1'] = '1 # the colomn number of VAR1'
    meta['CVAR2'] = '2 # the colomn number of VAR2'
    meta['NVAR1'] = 'NaN # NAXIS1'
    meta['NVAR2'] = '%d # NAXIS2'%(len(list_beta)*len(list_T_dust))
    meta['NPAR'] = '3 # Parameter 1 2 3 ...'
    meta['TPAR1'] = 'beta'
    meta['TPAR2'] = 'T_dust'
    meta['TPAR3'] = 'L_dust # integrated over 1-1000um. Final L_IR = TPAR3 * a * dL**2 / (1.+z).'
    meta['NPAR1'] = '%d # number of beta Values'%(len(list_beta))
    meta['NPAR2'] = '%d # number of T_dust values'%(len(list_T_dust))
    meta['NPAR3'] = '1 # non-independent'
    meta['CPAR1'] = '3'
    meta['CPAR2'] = '4'
    meta['CPAR3'] = '5'
    
    lambda_um = []
    Snu_Jy = []
    beta = []
    T_dust = []
    L_dust = []
    for i in range(len(list_beta)):
        for j in range(len(list_T_dust)):
            print('list_beta[%d/%d], list_T_dust[%d/%d]'%(i+1,len(list_beta),j+1,len(list_T_dust)))
            SED = get_SED_MBB(beta = list_beta[i], T_dust = list_T_dust[j], verbose = False)
            lambda_um.extend(SED['lambda_um'])
            Snu_Jy.extend(SED['Snu_Jy'])
            beta.extend(SED['Snu_Jy']*0.0 + list_beta[i])
            T_dust.extend(SED['Snu_Jy']*0.0 + list_T_dust[j])
            inte_val = 0.0
            for k in range(1,len(SED['lambda_um'])):
                if SED['lambda_um'][k-1] >= 8.0 and SED['lambda_um'][k] <= 1000.0:
                    freq_step = np.abs(2.99792458e5/(SED['lambda_um'][k])-2.99792458e5/(SED['lambda_um'][k-1])) # GHz
                    inte_step = (SED['Snu_Jy'][k]+SED['Snu_Jy'][k-1])/2.0*1e3 * freq_step / 40.31970 # 1 Lsun Mpc-2 = 40.31970 mJy GHz
                    inte_val += inte_step
            L_dust.extend(SED['Snu_Jy']*0.0 + inte_val)
    # 
    meta['NVAR1'] = '%d # NAXIS1'%(len(SED['lambda_um']))
    # 
    otb = Table()
    otb['lambda_um'] = lambda_um
    otb['Snu_Jy'] = Snu_Jy
    otb['beta'] = beta
    otb['T_dust'] = T_dust
    otb['L_dust'] = L_dust
    
    otb['lambda_um'].format = '%0.6e'
    otb['Snu_Jy'].format = '%0.6e'
    otb['beta'].format = '%0.1f'
    otb['T_dust'].format = '%0.1f'
    otb['L_dust'].format = '%0.6e'
    
    pprint(meta)
    #otb.meta = meta
    #with open('lib.MBB.SED.json', 'w') as fp:
    #    json.dump(meta, fp, indent=4, sort_keys=False)
    
    otb.write('lib.MBB.SED', format='ascii.fixed_width', delimiter=' ', overwrite=True)
    print('Output to "lib.MBB.SED"!')
























