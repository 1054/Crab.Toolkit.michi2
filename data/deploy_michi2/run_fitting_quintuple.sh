#!/bin/bash
#

echo "macro read rShift.sm redShift_maskout_LowSNR $1" | sm

./michi2_v03 -obs flux_restframe.dat \
             -lib FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED \
                  MullaneyAGN.Single.lib.SED \
                  DL07.HiExCom.SPAH.lib.SED \
                  DL07.LoExCom.SPAH.lib.SED \
                  RadioPowerlaw.Single.lib.SED \
             -out fit_quintuple.out \
             -constraint LIB3 INDEX EQ LIB4 INDEX
