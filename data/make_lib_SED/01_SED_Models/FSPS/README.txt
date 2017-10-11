http://people.ucsc.edu/~conroy/FSPS.html

Flexible Stellar Population Synthesis

FSPS is a flexible SPS package that allows the user to compute simple stellar populations (SSPs) for a range of IMFs and metallicities, and for a variety of assumptions regarding the morphology of the horizontal branch, the blue straggler population, the post--AGB phase, and the location in the HR diagram of the TP-AGB phase.

From these SSPs the user may then generate composite stellar populations (CSPs) for a variety of star formation histories (SFHs) and dust attenuation prescriptions.  As of v2.3 FSPS includes self-consistent incorporation of the Draine & Li 2007 dust emission spectra.  Outputs include the `observed' spectra and magnitudes of the SSPs and CSPs at arbitrary redshift.  In addition to these fortran routines, several IDL routines are provided that allow easy manipulation of the output.





$ cd /Users/dliu/Programming/sps/fsps/

$ svn checkout http://fsps.googlecode.com/svn/trunk/ fsps












FSPS was designed with the intention that the user would make full use of the provided fortran routines.  However, the full FSPS package is quite large, and requires some time for the user to become familiar with all of the options and syntax.  Some users may only need SSPs for a range of metallicities and IMFs.  For such users, standard SSP sets for several IMFs, evolutionary tracks, and spectral libraries are available here (http://people.ucsc.edu/~conroy/model_SSPs.html).  

* http://people.ucsc.edu/~conroy/model_SSPs.html

The filename specifies the IMF (Kroupa, Chabrier, or Salpeter), spectral library (BaSeL or Miles), and stellar evolutionary code (Padova).  Within each tarball you will find one file for each metallicity, where the metallicity is also part of the filename (in the units such that Z0.0190 is solar metallicity).

Each *.spec file contains a header followed by information for the spectral evolution of the SSP as a function of time.  There are two lines for each age.  The first line specifies the age in log(yrs), mass in log(Msun), bolometric luminosity in log(Lsun) and a dummy number.  The second line specifies the flux in fnu units.   See the manual for details regarding units.  

Within this directory you will also find the files basel.lambda and miles.lambda that list the wavelengths corresponding to the fluxes in the *.spec files.  The number of wavelength points in these files corresponds to the number of spectral elements within the *.spec files.  Wavelengths are in angstroms. 

Finally, this directory also contains an IDL routine (read_spec.pro) that will read the *.spec files into an IDL structure.



$ cd /Users/dliu/Programming/sps/fsps/

$ svn checkout http://fsps.googlecode.com/svn/trunk/SSP/ fsps_ssp



