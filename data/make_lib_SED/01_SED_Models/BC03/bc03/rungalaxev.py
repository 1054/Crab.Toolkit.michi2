#!/usr/bin/python
######################################################################
"""
Documentation
AMSR Nov 2007
additional tweaks pah

[] denotes a default
=[] denotes fixed
[/] denotes items in a restricted list.  The first option is the default
unless otherwise specfied (if not always).
[ , ] denotes a restricted range; default given separately
This version always computes its own input models.
All ages in Gyr

BASE NAME FOR OUTPUTS

FILENAME () [GalaxevOut]                        # Produces FILENAME plus
                                                   # .1ABmag .1color .2color .3color .4color .ised


INPUTS
Web page for Galaxev http://www2.iap.fr/users/charlot/bc2003/ (see documentation for references)

IMF [salpeter/chabrier]                            # Type of Initial Mass Function
RESOLUTION [lr/hr]                                 # Low/High resolution: up to 1221/6990 wavelength points 
TRACKS [Padova1994/Padova2000]                     # Stellar population evolutionary tracks
       if Padova1994
          METALLICITY [0.0001/0.0004/0.004/0.008/0.02/0.05] [0.02] # Default 0.02 is Solar metallicity
       if Padova2000
          METALLICITY [0.0004/0.001/0.004/0.008/0.019/0.03] [0.019] 
DUST [N/Y]                                          # Incude dust in models?
       if Y
          TAU_V [1.0]                               # Total effective attenuation optical depth
          MU    [0.3]                               # Fraction of tau_V arising from the ambient ISM
SFRTYPE [exponential/SSP/singleburst/constant/delayed/linear_decr]
       if SSP                                       # 0 Instantaneous burst of star formation
       if exponential                               # 1 SFR exponentially declining 
          TAU [0.01,2] [0.1]                        # e-folding timescale of decline TAU Gyr
                                                    # The code allows an alternative parameter MU_SFR
                                                    # but this is exactly determined by TAU, see README
                                                    # - let me know if this is wanted 
          GAS_RECYCLE [N/Y]                         # Recycle gas ejected by stars?
          if Y
             EPSILON [0.5]                          # Fraction of gas recycled
                                                    # (<0 implies wind, > 1 implies infall)
          TCUT  [20]                                # Set SFR back to 0 after TCUT Gyr     
       if singleburst                               # 2 Single burst of star formation of finite length
          TAU [0.05]                                # Duration in Gy
       if constant                                  # 3 Constant SFR 
          SFR [1]                                   # SFR in Msun/yr
          TCUT  [20]                                # Set SFR back to 0 after TCUT Gyr
       if delayed                                   # 4 SFR has delayed peak
          TAU [1]                                   # Time in Gy of max SFR in Msun/yr
          TCUT  [20]                                # Set SFR back to 0 after TCUT Gyr
       if linear_decr                               # 5 SFR is linearly decreasing
          TAU [1]                                   # SFR = 0 at time TAU Gyr
          TCUT  [20]                                # Set SFR back to 0 after TCUT Gyr
###################################################################################################             
W1 [91,1.600E+06] [100]                             # Lower limit of wavelength range (Angstrom)
W2 [91,1.600E+06] [1.2E+06]                         # Upper limit of wavelength range (Angstrom) 
UNITS [perAA, perHz]                                # Luminosity (Lsun) per spectral unit; covert to  perAA: W1=W1; perHz: W1=-1*W1
NORMALISATION [NONE/WAVE/Fz]                        # to allow the options to be presented if needed
       if WAVE                                      # Normalise spectrum to luminosity F0 at wavelength W0
          F0 [1]                                    
          W0 [W1, W2]
       if Fz                                        # Normalise spectrum such that F0 is the flux (Lsun/chosen spectral unit) from a source at redshift z in FILTER 
          F0 [1]
          z [1]
          FILTER [ 23 Johnson R ]                   # Filter No. N=23 is Johnson R; link to filters.log (omit 1st column?)
AGE  [0,20] [5,10,15,20]                            # 1-24 separate ages at which to present output (Gyr) **?? how separated?
####################################################################################################
COSMO_EVO [N/Y]                                     # Include cosmological evolution?
       if [Y]
          H0 [70]                                   # Hubble constant km/s/Mpc
          Omega [0.3]                               # Cosmological parameters
          Omega_lambda [0.7]                        # Cosmological parameters
          TG [13.0]                                 # Present age of galaxy Gyr
          FILTER1 [ 23 Johnson R ]                  # One or 2 filters, see FILTER above
          FILTER2 [109 2MASS J ]                    

"""
############################################################

import os
import sys


# Enter your parameters here

#input = ['./galaxevtest.py', 'FILENAME=GalaxevOut', 'IMF=salpeter', 'RESOLUTION=lr', 'TRACKS=Padova1994', 'METALLICITY=0.02', 'DUST=N', 'TAU_V=', 'MU=', 'SFRTYPE=exponential', 'TAU=0.1', 'GAS_RECYCLE=N', 'EPSILON=', 'SFR=1', 'TCUT=20', 'W1=100', 'W2=1.2E+06', 'UNITS=perAA', 'NORMALISATION=NONE', 'F0=', 'W0=', 'z=', 'FILTER= 23 Johnson R ', 'AGE=5,10,15,20', 'COSMO_EVO=N', 'H0=70', 'Omega=0.3', 'Omega_lambda=0.7', 'TG=13.0', 'FILTER1= 23 Johnson R ', 'FILTER2=109 2MASS J '] 
galaxevdir = os.path.dirname(sys.argv[0])+ os.sep
input = ['./galaxevtest.py', 'FILENAME=GalaxevOut', 'IMF=salpeter', 'RESOLUTION=lr', 'TRACKS=Padova1994', 'METALLICITY=0.02', 'DUST=N', 'TAU_V=3.0', 'MU=0.1', 'SFRTYPE=constant', 'TAU=0.1', 'GAS_RECYCLE=Y', 'EPSILON=0.5', 'SFR=1', 'TCUT=20', 'W1=100', 'W2=1.2E+06', 'UNITS=perHz', 'NORMALISATION=Fz', 'F0=1.0', 'W0=5000', 'z=1.0', 'FILTER= 23 Johnson R ', 'AGE=5,10,15,20', 'COSMO_EVO=Y', 'H0=75', 'Omega=0.31', 'Omega_lambda=0.69', 'TG=12.0', 'FILTER1= 23 Johnson R ', 'FILTER2=109 2MASS J ']

hashtemplate = {'FILENAME': 'GalaxevOut', 'IMF': 'salpeter', 'RESOLUTION': 'lr', 'TRACKS': 'Padova1994', 'METALLICITY': 0.02, 'DUST': 'N', 'TAU_V': 1.0, 'MU': 0.3, 'SFRTYPE': 'exponential', 'TAU': 0.1, 'GAS_RECYCLE': 'N', 'EPSILON': 0.5, 'SFR': 1.0, 'TCUT': 20.0,  'W1': 100.0, 'W2': 1.2E+06, 'UNITS': 'perAA', 'NORMALISATION': 'NONE', 'F0': 1.0, 'W0': 5000, 'z': 1.0, 'FILTER': ' 23 Johnson R ', 'AGE': '5,10,15,20','COSMO_EVO': 'N', 'H0': 70, 'Omega': 0.3, 'Omega_lambda': 0.7, 'TG': 13.0, 'FILTER1': ' 23 Johnson R ', 'FILTER2': '109 2MASS J ', 'REDSHIFT_EVOL': 'col'}

inhash = {}
for i in range(1, len(input)):
    x = input[i].split('=')
    inhash[x[0]] = x[1]

hashtemplate.update(inhash)
inhash = hashtemplate

# Set up output filenanes TBC

# Validate free-form inputs TO FOLLOW

# Keys for metallicity labels

P94 = {'0.0001': '_m22', '0.0004': '_m32', '0.004': '_m42', '0.008': '_m52', '0.02': '_m62', '0.05': '_m72'}
P00 = {'0.0004': '_m122', '0.001': '_m132', '0.004': '_m142', '0.008': '_m152', '0.019': '_m162', '0.03': '_m172'}

# Keys for SFR type

SFRType = {'exponential': 1, 'SSP': 0, 'singleburst': 2, 'constant': 3, 'delayed': 5, 'linear_decr': 5}

# Keys for filter codes

filtercodes = {" 13 Koo-Kron U+ ": 13, " 23 Koo-Kron J+ ": 23, " 25 Koo-Kron F+ ": 25, " 28 Koo-Kron N+ ": 28, " 15 Koo-Kron R ": 15, " 43 BJ (photographic) ": 43, " 23 RF (photographic) ": 23, " 13 Koo-Kron U+ ": 13, " 25 Koo-Kron J+ ": 25, " 25 Koo-Kron F+ ": 25, " 23 Koo-Kron N+ ": 23, " 24 Buser U ": 24, " 40 Buser B2 ": 40, " 40 Buser B3 ": 40, " 54 Buser V ": 54, " 13 M+S U ": 13, " 22 M+S B ": 22, " 25 M+S V ": 25, " 43 S+S B ": 43, " 38 S+S V ": 38, " 43 S+S R ": 43, " 33 ST-UV14 ": 33, " 42 ST-UV17 ": 42, " 53 ST-UV22 ": 53, " 64 ST-UV27 ": 64, " 33 OAO-UV1 ": 33, " 47 OAO-UV2 ": 47, " 59 OAO-UV3 ": 59, " 68 OAO-UV4 ": 68, " 58 OAO-UV5 ": 58, " 70 OAO-UV6 ": 70, " 23 Johnson R ": 23, " 26 Johnson I ": 26, " 31 Johnson J ": 31, " 17 Johnson K ": 17, " 14 Johnson L ": 14, " 51 Butcher r ": 51, " 36 Butcher i ": 36, " 31 Butcher-Oemler R ": 31, " 31 Butcher-Oemler R ": 31, " 23 Basel u ": 23, " 46 Basel g ": 46, " 22 Basel r ": 22, " 33 UKIRT H ": 33, " 53 R. S. Ellis U(PE) ": 53, "105 R. S. Ellis J ": 105, " 57 R. S. Ellis R ": 57, "111 R. S. Ellis N ": 111, "189 MacKay+Hall KG3 ": 189, "126 MacKay+Hall I ": 126, " 29 Gunn g Palomar ": 29, " 47 Gunn r Palomar ": 47, " 75 Gunn i Palomar ": 75, " 81 Gunn z Palomar ": 81, " 21 IR J Palomar ": 21, " 27 IR H Palomar ": 27, " 33 IR K Palomar ": 33, " 33 IR K Palomar ": 33, " 33 IR K Palomar ": 33, " 33 IR K Palomar ": 33, " 40 Tyson J ": 40, " 41 Tyson R ": 41, " 65 Tyson I ": 65, " 4 ANS 1550 Wide ": 4, " 18 ANS 1800 ": 18, " 4 ANS 2200 ": 4, " 4 ANS 2500 ": 4, " 23 ANS 3300 ": 23, " 13 c.U Lilly+Cowie ": 13, " 61 c.I Lilly+Cowie ": 61, " 35 IRAS 12 um ": 35, " 28 IRAS 25 um ": 28, " 46 IRAS 60 um ": 46, " 32 IRAS 100 um ": 32, " 20 H B+B ": 20, " 21 J B+B ": 21, " 23 K B+B ": 23, " 21 L 3.5 um B+B ": 21, " 20 L' 3.8 um B+B ": 20, " 19 M B+B ": 19, "200 KPNO 11.1um FWHM=1.72 ": 200, "200 KPNO 8.4 um FWHM=1.7 ": 200, " 7 Burstein+ 1250-1850A ": 7, " 65 Cousins R @ V=15. ": 65, " 39 Cousins I @ V=15. ": 39, " 65 W+C K' ": 65, " 71 W+C I ": 71, " 19 Washington C ": 19, " 18 Washington M ": 18, " 22 Washington T1 ": 22, " 18 Washington T2 ": 18, " 19 WFPC2 F300W ": 19, " 22 WFPC2 F450W ": 22, " 30 WFPC2 F606W ": 30, " 33 WFPC2 F814W ": 33, "327 WFPC2 F300W red leak ": 327, " 9 UV-1 2200-2400 A ": 9, " 6 UV-2 2640-2750 A ": 6, " 9 UV-3 2950-3150 A ": 9, "763 J-ESO ": 763, "557 H-ESO ": 557, "335 K-ESO ": 335, "160 L-ESO ": 160, "492 LW2-ISO 7 um ": 492, "311 I DeNIS ": 311, "251 J DeNIS ": 251, "300 Ks DeNIS ": 300, "400 CTIO Ks +NICMOS-3+atm ": 400, " 41 corr Cousins I @ V=1 ": 41, "275 u* MEGACAM (CFH) ": 275, "453 g' MEGACAM (CFH) ": 453, "503 r' MEGACAM (CFH) ": 503, "943 I' MEGACAM (CFH) ": 943, "745 z' MEGACAM (CFH) ": 745, " 47 SDSS u ": 47, " 89 SDSS g ": 89, " 75 SDSS r ": 75, " 89 SDSS i ": 89, "141 SDSS z ": 141, "109 2MASS J ": 109, " 58 2MASS H ": 58, " 78 2MASS Ks ": 78}

# Set environment etc. 
# setenv bc03 /scratch/hecate_1/amsr/ASTROGRID/SOFTWARE/GALAXEV/bc03/src
# cd $bc03
# source ./.bc_cshrc or equivalent - this contains:
## In your .cshrc or .login file, define bc03 as the directory which
## contains the GALAXEV programs (setenv bc03 /full_path_to_GALAXEV_directory)
##
#setenv FILTERS            $bc03/FILTERBIN.RES
#setenv A0VSED             $bc03/A0V_KURUCZ_92.SED
#setenv RF_COLORS_ARRAYS   $bc03/RF_COLORS.filters
#####################################################
#alias  add                'csh $bc03/add_bursts.sh'
#alias  csp                'csh $bc03/csp_galaxev.sh'
#alias  vdisp              'csh $bc03/vel_disp.sh'
#alias  cmev               'csh $bc03/cm_evolution.sh'
#alias  dgr                'csh $bc03/downgrade_resolution.sh'
#alias  gpl                '$bc03/galaxevpl'
#alias  zmag               '$bc03/zmag'
# set working directory
outdir = "./"

# Open shellscript for csp
cspshell = outdir + 'mycsp.sh'
mycsp = open(cspshell, 'w')
os.chmod(cspshell, 0777)
print >> mycsp, '#!/bin/csh \n'
print >> mycsp, 'setenv bc03 "'+galaxevdir+'src/"' 
print >> mycsp, 'setenv FILTERS $bc03/FILTERBIN.RES'

print >> mycsp, galaxevdir +'src/csp_galaxev.sh <<EOF'

# Set up path to tracks e.g. bc03/models/
# tracks are of form 
# Padova1994/salpeter/bc2003_hr_m32_salp_ssp.ised_ASCII

if (inhash['IMF'] == 'salpeter'):
    myimf = '_salp_ssp.ised_ASCII'

if (inhash['IMF'] == 'chabrier'):
    myimf = '_chab_ssp.ised_ASCII'

filetracksdir = galaxevdir + "models/"

if (inhash['TRACKS'] == 'Padova1994'):
    myfiletrack = filetracksdir + inhash['TRACKS'] + "/" + inhash['IMF'] + "/bc2003_" + inhash['RESOLUTION'] + P94[inhash['METALLICITY']] + myimf

if (inhash['TRACKS'] == 'Padova2000'):
    myfiletrack = filetracksdir + inhash['TRACKS'] + "/" + inhash['IMF'] + "/bc2003_" + inhash['RESOLUTION'] + P00[inhash['METALLICITY']] + myimf

localtrackfile = os.path.basename(myfiletrack)
os.symlink(myfiletrack, localtrackfile)
os.system(galaxevdir+"src/bin_ised " + localtrackfile )
print >> mycsp, localtrackfile

# calculate spectra for tracks

# set associated inputs

print >>  mycsp, inhash['DUST']
if (inhash['DUST'] == 'Y'):
    print >>  mycsp, inhash['TAU_V']
    print >>  mycsp, inhash['MU']
    
MYSFR = SFRType[inhash['SFRTYPE']]
print >> mycsp, MYSFR
    
if (MYSFR == 1):
    print >> mycsp, inhash['TAU']
    print >> mycsp, inhash['GAS_RECYCLE']
    if (inhash['GAS_RECYCLE'] == 'Y'):
        print >> mycsp, inhash['EPSILON']
    
    print >> mycsp, inhash['TCUT']

if (MYSFR == 2):
    print >> mycsp, inhash['TAU']

if (MYSFR == 3):
    print >> mycsp, inhash['SFR']
    print >> mycsp, inhash['TCUT']

if (MYSFR == 4):    
    print >> mycsp, inhash['TAU']
    print >> mycsp, inhash['TCUT']

if (MYSFR == 5):
    print >> mycsp, inhash['TAU']
    print >> mycsp, inhash['TCUT']

print >> mycsp, inhash["FILENAME"]
print >> mycsp, 'EOF' 
mycsp.close()

# run - check path****

os.system("./mycsp.sh")

####
# Open shellscript for gpl
gplshell = outdir + 'mygpl.sh'
mygpl = open(gplshell, 'w')
os.chmod(gplshell, 0777)
print >> mygpl, '#!/bin/csh'
print >> mygpl, 'setenv bc03 "'+galaxevdir+'src/"' 
print >> mygpl, 'setenv FILTERS $bc03/FILTERBIN.RES'

print >> mygpl, galaxevdir +'src/galaxevpl <<EOF'

print >> mygpl, inhash["FILENAME"]

# set up extraction parameters
W = float(inhash['W1'])
if (inhash['UNITS'] == 'perHz'):
    W = -1.0*W

W = str(W) + ',' + str(inhash['W2']) 

if (inhash['NORMALISATION'] == 'WAVE'):
    W = W + ',' + str(inhash['W0']) + ',' + str(inhash['F0'])

if (inhash['NORMALISATION'] == 'Fz'):
    N = ',' + str(-1.0 * float(filtercodes[inhash['FILTER']])) + ','
    W = W + N + ',' + str(inhash['F0']) + ',' + str(inhash['z'])

print >> mygpl, W    
print >> mygpl, inhash['AGE']
print >> mygpl, "\n"
print >> mygpl, 'EOF'
mygpl.close()

# run - check path****

os.system("./mygpl.sh")

######
# Open shellscript for cmev

if (inhash['COSMO_EVO'] == 'Y'):
    cmevshell = outdir + 'mycmev.sh'
    mycmev = open(cmevshell, 'w')
    os.chmod(cmevshell, 0777)
    print >> mycmev, '#!/bin/csh \n'
    print >> mycmev, 'setenv bc03 "'+galaxevdir+'src/"' 
    print >> mycmev, 'setenv FILTERS $bc03/FILTERBIN.RES'

    print >> mycmev, galaxevdir +'src/cm_evolution.sh <<EOF'
    print >> mycmev,  str(inhash['H0']) + ',' + str(inhash['Omega']) + ',' + str(inhash['Omega_lambda'])
    print >> mycmev, inhash['TG']
    print >> mycmev, inhash["FILENAME"]
    print >> mycmev, filtercodes[inhash['FILTER1']], ',', filtercodes[inhash['FILTER2']]
    print >> mycmev, '\n'
    print >> mycmev, 'EOF'
    mycmev.close()

    # run - check path****

    os.system("./mycmev.sh")

#do some file renaming to make it easier to pick up the files
os.rename(inhash["FILENAME"]+".magnitude_F"+"%03d" % int(filtercodes[inhash["FILTER1"]]), inhash["FILENAME"]+".magnitude_FILTER1")
if inhash["REDSHIFT_EVOL"] == "col":
    os.rename(inhash["FILENAME"]+".magnitude_F"+"%03d" % int(filtercodes[inhash["FILTER2"]]), inhash["FILENAME"]+".magnitude_FILTER2")

