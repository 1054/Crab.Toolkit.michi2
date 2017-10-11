#!/bin/bash
# 

#ls outputs/FSPS.Padova.BaSeL.Z0.0190.lib.SED

echo "macro read make_Calzetti2000law.sm make_Calzetti2000law_BC03" | sm

zip lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.zip lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED.readme

