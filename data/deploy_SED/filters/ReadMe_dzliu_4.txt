# Downloading filter curves from https://cosmos.astro.caltech.edu/page/filterset
# then converting wavelength to microns.

sm << EOF


data filter_curve_UVISTA_Y.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UVISTA_Y.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_NEWFIRM_H1.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_KPNO_NEWFIRM_H1.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_NEWFIRM_H2.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_KPNO_NEWFIRM_H2.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_NEWFIRM_J1.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_KPNO_NEWFIRM_J1.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_NEWFIRM_J2.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_KPNO_NEWFIRM_J2.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_NEWFIRM_J3.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_KPNO_NEWFIRM_J3.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_NEWFIRM_Ks.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_KPNO_NEWFIRM_Ks.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_UVISTA_H.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UVISTA_H.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_UVISTA_J.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UVISTA_J.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_UVISTA_Ks.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UVISTA_Ks.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

quit
EOF

