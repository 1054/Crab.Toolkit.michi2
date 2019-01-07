#!/bin/bash
# 

set -e



cp "$HOME/Softwares/BC03/bc03_Updated_version_2016/bc03/out/out_constant_SFH_wave_flux_at_"*".txt" \
   "$HOME/Softwares/BC03/bc03_Updated_version_2016/bc03/out/out_constant_SFH_age_and_mass_tot.txt" \
    ./

./makelibSED_BC03.sh



./make_Calzetti2000law.sh



old_file_name="lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED"
new_file_name="lib.BC03.Padova1994.BaSeL.Z0.0190.ConstSFH.MultiAge.EBV.SED"
mv "$old_file_name" "$new_file_name"

if [[ $(uname) == Darwin ]]; then
    greadlink -f "$new_file_name" > "$new_file_name".readme
else
    readlink -f "$new_file_name" > "$new_file_name".readme
fi

zip "$new_file_name.zip" \
    "$new_file_name" \
    "$new_file_name".readme

