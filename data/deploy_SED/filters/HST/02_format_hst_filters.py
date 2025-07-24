#!/usr/bin/env python
#
import glob
import numpy as np

for filename in list(sorted(glob.glob('./HST_*.dat'))):
    outfilename = filename.replace('.dat', '.formatted.txt')
    table = np.loadtxt(filename, dtype={'names': ['wave_AA', 'transmission'], 'formats': ['f4', 'f4']})
    wave_um = np.linspace(np.min(table['wave_AA']), np.max(table['wave_AA']), 200) / 1e4
    transmission = np.interp(wave_um * 1e4, table['wave_AA'], table['transmission'])
    transmission[0] = 0.0
    transmission[-1] = 0.0
    with open(outfilename, 'w') as fp:
        fp.write('# microns transmission\n')
        for i in range(len(wave_um)):
            fp.write('{:12.5f} {:12.5e}\n'.format(wave_um[i], transmission[i]))

