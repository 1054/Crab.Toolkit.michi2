#!/bin/bash
#
sm << EOF
define No_Recompute_Radio 1
macro read pChisq.sm plotChisq \
"flux_obsframe.dat" \
"FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED" \
"DL07.HiExCom.SPAH.lib.SED" \
"DL07.LoExCom.SPAH.lib.SED" \
"RadioPowerlaw.Single.lib.SED" \
"fit_quartered.out"
quit
EOF

ln -fs fit_quartered.out fit.out
ln -fs run_plotting_quartered.sm run_plotting.sm

sm << EOF
define No_Recompute_Radio 1
macro read rChisq.sm
calc_chisq_probability
calc_uncertainty_of_U
calc_uncertainty_of_LIR
calc_uncertainty_of_Mdust
calc_uncertainty_of_Mstar
EOF
