#!/bin/bash
#

if [[ $# -eq 0 ]]; then
    echo "Please input redshift!"; exit
fi

#echo "macro read rShift.sm redShift_maskout_LowSNR $1" | sm

./michi2_v04 -obs flux_obsframe.dat \
             -redshift $1 \
             -lib lib.BC03.Padova1994.BaSeL.Z0.0190.EBV.SED \
                  MullaneyAGN.Single.lib.SED \
                  lib.DL07.HiExCom.SED \
                  lib.DL07.LoExCom.SED \
                  RadioPowerlaw.Single.lib.SED \
             -out fit_5.out \
             -parallel 2 \
             -constraint LIB3 INDEX EQ LIB4 INDEX \
             -constraint LIB3 PAR2 GT VALUE 0.400000 \
             -constraint LIB3 PAR2 LT VALUE 0.500000 \
             -constraint LIB4 PAR2 GT VALUE 0.400000 \
             -constraint LIB4 PAR2 LT VALUE 0.500000 \
             -constraint LIB5 NORM EQ SED "vLv(8,1000)*1.061619121e-06" -filter filter.list


