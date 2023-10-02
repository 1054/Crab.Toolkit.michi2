# Downloading filter curves from https://cosmos.astro.caltech.edu/page/filterset
# then converting wavelength to microns.

sm << EOF

data filter_curve_2MASS_J_2mass.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_2MASS_J.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_2MASS_H_2mass.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_2MASS_H.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_2MASS_Ks_2mass.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_2MASS_Ks.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_CFHT_u_megaprime_sagem.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_CFHT_MegaPrime_u.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_CFHT_g_megaprime_sagem.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_CFHT_MegaPrime_g.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_CFHT_r_megaprime_sagem.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_CFHT_MegaPrime_r.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_CFHT_i_megaprime_sagem.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_CFHT_MegaPrime_i.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_CFHT_z_megaprime_sagem.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_CFHT_MegaPrime_z.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_CFHT_wircam_H.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_CFHT_wircam_H.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_CFHT_wircam_Ks.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_CFHT_wircam_Ks.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data filter_curve_KPNO_flamingos_Ks.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_KPNO_flamingos_Ks.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

quit
EOF

