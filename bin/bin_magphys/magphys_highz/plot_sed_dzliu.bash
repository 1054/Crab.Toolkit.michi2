#!/bin/bash
# 

export IDL_PATH="$IDL_PATH:$HOME/Cloud/Github/Crab.Toolkit.michi2/bin/bin_magphys/magphys_highz"
export USER_OBS="magphys_input_fluxes.dat"
export USER_FILTERS="magphys_input_filters.dat"
export magphys="$HOME/Softwares/magphys/magphys_highz"

idl -e plot_sed -args "best-fit_SED"

ps2pdf -dEPSCrop "best-fit_SED.ps" "best-fit_SED.pdf"
