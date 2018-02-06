
# 
# 2018-02-06
# 



# 
# Compile FSPS
# 
git clone https://github.com/cconroy20/fsps
cd src
vim Makefile 
#-> change "F90 = gfortran" to "F90 = gfortran-mp-5"
vim sps_vars.f90
#-> change "#define MILES 1" to "#define MILES 0"
#-> change "#define BASEL 0" to "#define BASEL 1"
#-> change "REAL(SP) :: om0=0.27, ol0=0.73, H0=72." to "REAL(SP) :: om0=0.27, ol0=0.73, H0=73."
#-> change to "INTEGER :: add_dust_emission=0"
#-> change to "INTEGER :: add_agn_dust=0"
#-> change to "INTEGER :: add_agb_dust_model=0"
#-> change to "INTEGER :: dust_type=2"
vim autosps.f90
make clean
make



# 
# Run FSPS with 
#               Chabrier IMF
#               tau-declining SFH with no single busrt
#               1.0 solar metallicity
# 
export PATH="$HOME/Cloud/Github/fsps/src:$PATH"
export SPS_HOME="$HOME/Cloud/Github/fsps"
autosps.exe
######
 enter IMF [0-5; def:0]:
  (0=Salpeter, 1=Chabrier 2003, 2=Kroupa 2001, 3=van Dokkum 2008, 4=Dave 2008, 5=tabulated)
1
 ---> Using IMF 1
 Specify SFH [0-2, def:0]
 (0=SSP, 1=CSP, 2=tabulated)
1
 ---> Computing a CSP
 input parameters for CSP: tau, const, age, fburst, tburst
  - tau in Gyr
  - const as fraction of mass formed in constant component
  - age of the system in Gyr.  i.e. results span the time 0<t<age
  - fburst as fraction of mass formed in an instantaneous burst
  - tburst as time of burst, with tburst<age
1.0
0.0
0.0
0.0
99.0
 ---> (tau const age fburst tburst)=(  1.00  0.00  0.00  0.00 99.00)
 enter metallicity [1-22; def:20]:
20
 Include default dust model? [yes/no, def:no]
 (default: tau1=1.0, tau2=0.3, MW extinction)
no
 ---> tau1= 0.00, tau2= 0.00
 Enter filename [def: "CSP.out"]
 ---> Output filename: CSP.out                                                                                             
 ---> Running model.......
######
cd $SPS_HOME/OUTPUTS/
cat CSP.out.spec | head -n 10 | tail -n 2 > output_row_w.txt
cat CSP.out.spec | tail -n +11 > output_row_f.txt
######
export IDL_PATH="$IDL_PATH:$(dirname $(pwd))/pro"
idl
res = read_spec('CSP.out.spec')
print, tag_names(res) ; AGEGYR LOGMASS LOGLBOL LOGSFR SPEC LAMBDA
print, size(res.LAMBDA)
print, size(res.SPEC)
log_ages = ALOG10(res.AGEGYR*1e9)
log_masses = res.LOGMASS
log_Lbols = res.LOGLBOL
log_SFRs = res.LOGSFR
n_spec = (size(res.LAMBDA,/DIM))[1]
n_wave = (size(res.LAMBDA,/DIM))[0]
flux = dblarr(n_spec+1,n_wave)
flux[0,*] = (res.LAMBDA)[*,0]
for i = 1, n_spec do flux[i,*] = (res.SPEC)[*,i-1]
CrabTablePrintF, "output_wave_fluxes.dat", flux
CrabTablePrintC, "output_spec_params.dat", log_ages, log_masses, log_Lbols, log_SFRs
exit





