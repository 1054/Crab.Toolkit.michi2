#!/bin/bash
# 

cd outputs/

echo "macro read \"../make_Calzetti2000law.sm\" make_Calzetti2000law_BC03" | sm

greadlink -f "lib.FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.EBV.SED" > "lib.FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.EBV.SED.readme"

zip "lib.FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.EBV.SED.zip" "lib.FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.EBV.SED" "lib.FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.EBV.SED.readme"

cd ../

