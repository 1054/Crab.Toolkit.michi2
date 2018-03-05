#!/usr/bin/env python
# 

import os
import sys
import numpy
import astropy
import astropy.io.ascii as asciitable
from copy import copy

Script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(Script_dir)


from michi2_read_lib_SEDs import *

table_stellar = asciitable.read('test_1/test_stellar.txt')
w_stellar = table_stellar['Wavelength_um'].data
f_stellar = table_stellar['Flux_density_mJy'].data
f_stellar_unattenuated = table_stellar['Flux_density_unattenuated_mJy'].data

z = 1.5

print('integrate_vLv stellar SED attenuated = %0.6e [Lsolar]'%(integrate_vLv(w_stellar, f_stellar, z)))
print('integrate_vLv stellar SED unattenuated = %0.6e [Lsolar]'%(integrate_vLv(w_stellar, f_stellar_unattenuated, z)))


