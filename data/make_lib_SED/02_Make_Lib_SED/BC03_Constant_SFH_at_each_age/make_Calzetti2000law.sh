#!/bin/bash
# 

echo "macro read make_Calzetti2000law.sm make_Calzetti2000law_BC03" | sm

lib_file_name="lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED"

# rm "lib.BC03.Padova1994.BaSeL.Z0.0190.ConstSFH.Age400Myr.SED"

if [[ $(uname) == Darwin ]]; then
    greadlink -f "$lib_file_name" > "$lib_file_name".readme
else
    readlink -f "$lib_file_name" > "$lib_file_name".readme
fi

#zip "$lib_file_name".zip \
#    "$lib_file_name" \
#    "$lib_file_name".readme

