#!/bin/bash
#

if [[ $# -eq 0 ]]; then
    echo "Please input redshift!"; exit
fi

echo "macro read rShift.sm redShift_maskout_LowSNR $1" | sm

./michi2_v03 -obs flux_restframe.dat \
             -lib FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED \
                  DL07.HiExCom.SPAH.lib.SED \
                  DL07.LoExCom.SPAH.lib.SED \
                  RadioPowerlaw.Single.lib.SED \
             -out fit_quartered.out \
             -constraint LIB2 INDEX EQ LIB3 INDEX
