wget https://raw.githubusercontent.com/GabrielaCR/AGNfitter/master/models/FILTERS/SPITZER/irac_ch1.res
wget https://raw.githubusercontent.com/GabrielaCR/AGNfitter/master/models/FILTERS/SPITZER/irac_ch2.res
wget https://raw.githubusercontent.com/GabrielaCR/AGNfitter/master/models/FILTERS/SPITZER/irac_ch3.res
wget https://raw.githubusercontent.com/GabrielaCR/AGNfitter/master/models/FILTERS/SPITZER/irac_ch4.res

sm

data irac_ch1.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_irac_ch1.dat" '%12.4e %15.6e\n' {wavelength_um transmission}

data irac_ch2.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_irac_ch2.dat" '%12.4e %15.6e\n' {wavelength_um transmission}

data irac_ch3.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_irac_ch3.dat" '%12.4e %15.6e\n' {wavelength_um transmission}

data irac_ch4.res read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_irac_ch4.dat" '%12.4e %15.6e\n' {wavelength_um transmission}

quit

