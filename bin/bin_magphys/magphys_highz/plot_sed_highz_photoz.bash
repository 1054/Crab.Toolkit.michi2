#!/bin/bash
# 
# Usage:
#    please change directory to where the magphys fitting results *.sed and *.fit are stored, 
#    then run this code 
# Example:
#    ssh -Y aida40110
#    cd /disk1/dzliu/Works/AlmaCosmos/Samples/20180720/Multi-wavelength_SEDs/ID_433236/SED_fitting_magphys_photoz/magphys_fitting/fit_1_with_datatable_photometry_magphys/
#    ls *.sed *.fit
#    ~/Softwares/magphys/magphys_photoz/plot_sed_highz_photoz.bash
#    gv 1.ps
#    evince 1_crop.pdf
# 

# readlinkfull
readlinkfull() {
    if [[ $# -gt 1 ]]; then if [[ "$1" == "-f" ]]; then shift; fi; fi
    DIR=$(echo "${1%/*}"); (cd "$DIR" && echo "$(pwd)/$(basename ${1})")
}

export magphys=$(readlinkfull $(dirname ${BASH_SOURCE[0]}))
export IDL_PATH="$IDL_PATH:+$magphys"
export FILTERS="$magphys/FILTERBIN.RES"
export OPTILIB="$magphys/OptiLIB_bc03_highz.bin"
export OPTILIBIS="$magphys/OptiLIBis_bc03_highz.bin"
export IRLIB="$magphys/InfraredLIB_highz.bin"
export USER_FILTERS="magphys_input_filters.dat"
export USER_OBS="magphys_input_fluxes.dat"

if [[ $# -gt 0 ]]; then 
idl -e plot_sed_highz_photoz -args $@
else
idl -e plot_sed_highz_photoz -args 1
ps2pdf -dEPSCrop 1.ps 1.pdf
pdfcrop -margin 10 1.pdf 1_crop.pdf
fi

