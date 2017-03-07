#!/bin/bash
#
sm << EOF
define No_Recompute_Radio 1
macro read pChisq.sm plotChisq \
"flux_obsframe.dat" \
"FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED" \
"MullaneyAGN.Single.lib.SED" \
"DL07.HiExCom.SPAH.lib.SED" \
"DL07.LoExCom.SPAH.lib.SED" \
"RadioPowerlaw.Single.lib.SED" \
"fit_quintuple.out"
quit
EOF

