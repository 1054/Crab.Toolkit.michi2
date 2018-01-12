#!/bin/bash
#
sm << EOF
macro read pChisq.sm plotChisq \
"flux_obsframe.dat" \
FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED \
DL07.LoExCom.SPAH.lib.SED \
RadioPowerlaw.Single.lib.SED \
fit.out
quit
EOF

#sm << EOF
#macro read pChisq.sm plotChisq \
#$(cat fit.out.info | grep "^OBS = "   | sed -e 's/.* *= *//g') \
#$(cat fit.out.info | grep "^LIB[0-9]" | sed -e 's/.* *= *//g') \
#$(cat fit.out.info | grep "^OUT = "   | sed -e 's/.* *= *//g')
#quit
#EOF

