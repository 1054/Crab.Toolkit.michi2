#!/bin/bash
# 
#export magphys="/Users/dzliu/Softwares/magphys/magphys_highz"
#export FILTERS="/Users/dzliu/Softwares/magphys/magphys_highz/FILTERBIN.RES"
#export OPTILIB="/Users/dzliu/Softwares/magphys/magphys_highz/OptiLIB_bc03_highz.bin"
#export OPTILIBIS="/Users/dzliu/Softwares/magphys/magphys_highz/OptiLIBis_bc03_highz.bin"
#export IRLIB="/Users/dzliu/Softwares/magphys/magphys_highz/InfraredLIB_highz.bin"
#export USER_FILTERS="/Users/dzliu/Work/DeepFields/Works_cooperated/2018_John_Silverman/SED_fitting_magphys_20180110/PACS-830/magphys_fitting/fit_1_with_flux_obsframe/magphys_input_filters.dat"
#export USER_OBS="/Users/dzliu/Work/DeepFields/Works_cooperated/2018_John_Silverman/SED_fitting_magphys_20180110/PACS-830/magphys_fitting/fit_1_with_flux_obsframe/magphys_input_fluxes.dat"
#export magphys_outdir="/Users/dzliu/Work/DeepFields/Works_cooperated/2018_John_Silverman/SED_fitting_magphys_20180110/PACS-830/magphys_fitting/fit_1_with_flux_obsframe"
#export michi2_outdir="../obj_1"

export magphys="/Users/dzliu/Softwares/magphys/magphys_highz"
export USER_FILTERS="magphys_input_filters.dat"
export magphys_outdir="."
export michi2_outdir="."


if [[ ! -f "SED_SUM_REDSHIFTED" ]]; then
sm << EOF
data "$michi2_outdir/SED_SUM"
read {w_um 1 f_mJy 2}
data "$michi2_outdir/redshift.txt"
read z 1
set w_um = w_um * (1+z)
print "SED_SUM_REDSHIFTED" '%15.6e %15.6e\n' {w_um f_mJy}
EOF
fi


if [[ ! -f 1.fit ]]; then
    cp "$magphys_outdir/1.fit" .
    cp "$magphys_outdir/1.sed" .
fi


idl -e plot_sed_overlaid_another_sed -args 1 "SED_SUM_REDSHIFTED"

ps2pdf -dEPSCrop 1.ps 1.pdf

pdfcrop -margin 5 1.pdf 1.pdf


