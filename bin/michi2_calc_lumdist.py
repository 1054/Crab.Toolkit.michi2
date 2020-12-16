#!/usr/bin/env python
# 

import os
import sys
import numpy
import shutil
import astropy

#from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)




####################################
#               MAIN               #
####################################

if len(sys.argv) <= 1:
    print('Usage: michi2_calc_lumdist.py redshift')
    sys.exit()


print('# %10s %10s'%('dL_Mpc', 'z'))

for i in range(1, len(sys.argv)):
    
    Redshift = float(sys.argv[i])
    
    dL = cosmo.luminosity_distance(Redshift).value
    
    print('  %10.2f %10.4f'%(dL, Redshift))

