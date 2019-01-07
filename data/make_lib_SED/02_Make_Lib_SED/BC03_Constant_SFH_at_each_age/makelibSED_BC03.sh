#!/bin/bash
# 

if [[ ! -f "./datatable_wave_AA_flux_Lsun_per_AA.txt" ]]; then
    echo "Error! \"./datatable_wave_AA_flux_Lsun_per_AA.txt\" was not found! Please copy it from BC03 outputs following ../ReadMe* files!"
    exit
fi


Age_value=200 # default 200Myr
if [[ $# -ge 1 ]]; then
    Age_value=$(echo "$1" | sed -e 's/Myr//g')
fi


idl << EOF
.r makelibSED_BC03
makelibSED_BC03, ${Age_value}
EOF


