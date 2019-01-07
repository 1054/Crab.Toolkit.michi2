#!/bin/bash
# 

#set -e


#Age_list=("200Myr" "300Myr" "400Myr" "500Myr" "600Myr" "700Myr" "800Myr" "900Myr" "1Gyr")
Age_list=("200Myr" "400Myr")

for Age_str in ${Age_list[@]}; do
    
    cp "$HOME/Softwares/BC03/bc03_Updated_version_2016/bc03/out/out_constant_SFH_wave_flux_at_${Age_str}.txt" \
        "./datatable_wave_AA_flux_Lsun_per_AA.txt"
    
    ./makelibSED_BC03.sh "$Age_str"
    
    
    
    ./make_Calzetti2000law.sh
    
    
    
    old_file_name="lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED"
    new_file_name="lib.BC03.Padova1994.BaSeL.Z0.0190.ConstSFH.Age${Age_str}.EBV.SED"
    rm "$old_file_name".readme
    mv "$old_file_name" "$new_file_name"
    
    if [[ $(uname) == Darwin ]]; then
        greadlink -f "$new_file_name" > "$new_file_name".readme
    else
        readlink -f "$new_file_name" > "$new_file_name".readme
    fi
    
done

