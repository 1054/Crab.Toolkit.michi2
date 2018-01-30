#!/bin/bash
# 

echo "macro read make_Calzetti2000law.sm make_Calzetti2000law_BC03" | sm

if [[ $(uname) == Darwin ]]; then
    greadlink -f lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED > lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.readme
else
    readlink -f lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED > lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.readme
fi

zip lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.zip lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.readme

