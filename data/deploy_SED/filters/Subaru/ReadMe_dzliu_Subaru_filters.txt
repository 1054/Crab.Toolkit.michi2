# Downloading filter curves from https://cosmos.astro.caltech.edu/page/filterset
# then converting wavelength to microns.

sm << EOF

data g_HSC.txt read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_HSC_g.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data i_HSC.txt read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_HSC_i.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data r_HSC.txt read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_HSC_r.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data y_HSC.txt read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_HSC_y.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data z_HSC.txt read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_HSC_z.dat" '%16.4e %15.6e\n' {wavelength_um transmission}



data B_subaru.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_B.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data V_subaru.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_V.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data g_subaru.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_g.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data r_subaru.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_r.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data i_subaru.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_i.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data z_subaru.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_z.dat" '%16.4e %15.6e\n' {wavelength_um transmission}





data IA427.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA427.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA464.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA464.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA484.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA484.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA505.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA505.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA527.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA527.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA574.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA574.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA624.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA624.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA679.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA679.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA709.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA709.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA738.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA738.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA767.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA767.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data IA827.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_IA827.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data NB711.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_NB711.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data NB816.SuprimeCam.pb read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_Subaru_SuprimeCam_NB816.dat" '%16.4e %15.6e\n' {wavelength_um transmission}


quit
EOF
