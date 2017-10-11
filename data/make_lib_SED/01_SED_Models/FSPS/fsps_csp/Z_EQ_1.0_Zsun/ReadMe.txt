

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
vim autosps.f90
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
0.0
 ---> (tau const age fburst tburst)=( 1.00 0.00 0.00 0.00 0.00)
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
idl
openr, lun_w, 'output_row_w.txt', /get_lun
point_lun, lun_w, 0
readf, lun_w, format = '(i0,i0)', n_spec, n_wave & f_spec = dblarr(n_wave)
readf, lun_w, f_spec & flux = dblarr(n_spec+1,n_wave) & flux[0,*] = f_spec
free_lun, lun_w
log_ages = [] & log_masses = [] & log_Lbols = [] & log_SFRs = [] & spec_id = []
log_age = 0.0 & log_mass = 0.0 & log_Lbol = 0.0 & log_SFR = 0.0
openr, lun, 'output_row_f.txt', /get_lun
for i = 1, n_spec do begin & $
readf, lun, log_age, log_mass, log_Lbol, log_SFR & $
readf, lun, f_spec & flux[i,*] = f_spec & $
log_ages = [log_ages, log_age] & $
log_masses= [log_masses, log_mass] & $
log_Lbols = [log_Lbols, log_Lbol] & $
log_SFRs = [log_SFRs, log_SFR] & $
endfor
free_lun, lun
CrabTablePrintF, "output_wave_fluxes.dat", flux
CrabTablePrintC, "output_spec_params.dat", log_ages, log_masses, log_Lbols, log_SFRs
exit





