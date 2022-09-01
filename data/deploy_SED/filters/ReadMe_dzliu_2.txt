# Downloading filter curves from http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
# then converting wavelength to microns.

# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=SLOAN/SDSS.uprime_filter
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=SLOAN/SDSS.gprime_filter
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=SLOAN/SDSS.rprime_filter
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=SLOAN/SDSS.iprime_filter
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=SLOAN/SDSS.zprime_filter
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=UKIRT/UKIDSS.Z
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=UKIRT/UKIDSS.Y
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=UKIRT/UKIDSS.J
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=UKIRT/UKIDSS.H
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=UKIRT/UKIDSS.K
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=WISE/WISE.W1
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=WISE/WISE.W2
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=WISE/WISE.W3
# wget http://svo2.cab.inta-csic.es/theory/fps3/getdata.php?format=ascii&id=WISE/WISE.W4

sm << EOF

data SLOAN_SDSS.uprime_filter.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_SDSS_uprime.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data SLOAN_SDSS.gprime_filter.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_SDSS_gprime.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data SLOAN_SDSS.rprime_filter.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_SDSS_rprime.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data SLOAN_SDSS.iprime_filter.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_SDSS_iprime.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data SLOAN_SDSS.zprime_filter.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_SDSS_zprime.dat" '%16.4e %15.6e\n' {wavelength_um transmission}


data UKIRT_UKIDSS.Z.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UKIDSS_Z.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data UKIRT_UKIDSS.Y.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UKIDSS_Y.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data UKIRT_UKIDSS.J.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UKIDSS_J.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data UKIRT_UKIDSS.H.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UKIDSS_H.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data UKIRT_UKIDSS.K.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_UKIDSS_K.dat" '%16.4e %15.6e\n' {wavelength_um transmission}


data WISE_WISE.W1.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_WISE1.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data WISE_WISE.W2.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_WISE2.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data WISE_WISE.W3.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_WISE3.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

data WISE_WISE.W4.dat read {wavelength_A 1.f transmission 2.f}
set wavelength_um = wavelength_A/1e4
print "filter_curve_WISE4.dat" '%16.4e %15.6e\n' {wavelength_um transmission}

quit
EOF
