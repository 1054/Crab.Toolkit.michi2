#!/bin/bash
# 

for Age_value in "200Myr"; do
#  "300Myr" "400Myr" "500Myr" "600Myr" "700Myr" "800Myr" "900Myr" "1Gyr"

cp  "../../01_SED_Models/BC03/bc03/out/out_constant_SFH_wave_flux_at_${Age_value}.txt" \
    'datatable_wave_AA_flux_Lsun_per_AA.txt'

idl -e 'makelibSED_BC03'

cp "lib.BC03.Padova1994.BaSeL.Z0.0190.SED" "lib.BC03.Padova1994.BaSeL.Z0.0190.AGE${Age_value}.SED"

done

