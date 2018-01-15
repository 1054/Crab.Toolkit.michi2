#!/bin/bash
# 

echo "macro read make_Calzetti2000law.sm make_Calzetti2000law_BC03" | sm

greadlink -f "lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED" > "lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.readme"

zip lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.zip "lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED" "lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.readme"

