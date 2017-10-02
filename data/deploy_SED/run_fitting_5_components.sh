#!/bin/bash
#

if [[ $# -eq 0 ]]; then
    echo "Please input redshift!"; exit
fi

#echo "macro read rShift.sm redShift_maskout_LowSNR $1" | sm

./michi2_v04 -obs flux_obsframe.dat \
             -redshift $1 \
             -lib FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED \
                  MullaneyAGN.Single.lib.SED \
                  DL07.HiExCom.lib.SED \
                  DL07.LoExCom.lib.SED \
                  RadioPowerlaw.Single.lib.SED \
             -out fit_5.out \
             -constraint LIB3 INDEX EQ LIB4 INDEX \
             -constraint LIB5 NORM EQ SED "vLv(8,1000)*1.061619121e-06" -filter filter.list


