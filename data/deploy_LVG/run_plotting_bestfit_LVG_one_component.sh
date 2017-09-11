#!/bin/bash
#
sm << EOF
macro read pChisq_LVG.sm plotChisq \
"flux_co.dat" \
"lib_z_1.500_co_dvddr_5.0_dVDoW_50.0_Faster.lvg" \
"fit_one_component.out"
quit
EOF

