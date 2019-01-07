#!/bin/bash
# 

#~/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_filter_flux_2sigma_fit_infrared_upper_limits.py \
#"ID1929_z5.85.txt" \
#"datatable_for_SED_fitting.txt"

cp ~/Cloud/Github/Crab.Toolkit.michi2/data/make_lib_SED/02_Make_Lib_SED/ModifiedBlackbody/lib.dust.MBB.SED \
.


cp ~/Cloud/Github/Crab.Toolkit.michi2/data/make_lib_SED/02_Make_Lib_SED/Radio_mm_FIR_Lines/lib.Radio.mm.FIR.lines.SED \
.


rm -rf "fit_with_CO.out" \
"fit_with_CO.out.info" \
"results_fit_with_CO/"


~/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_v05/michi2_v05_mac \
-obs "datatable_for_SED_fitting.txt" \
-lib "lib.dust.MBB.SED" \
     "lib.Radio.mm.FIR.lines.SED" \
-constraint "LIB2_NORM = LIB1_NORM * LIB1_PAR3 * 55663.9^2 / (1.+5.85) / LIB2_PAR5 * LIB2_PAR4^2 / 55663.9^2 * (1.+5.85)" \
-filter "filters/Herschel/PACS_70.dat" \
        "filters/Herschel/PACS_100.dat" \
        "filters/Herschel/PACS_160.dat" \
        "filters/Herschel/SPIRE_250.dat" \
        "filters/Herschel/SPIRE_350.dat" \
        "filters/Herschel/SPIRE_500.dat" \
-redshift 5.85 \
-sampling 1000 \
-out "fit_with_CO.out" \
-parallel 2
# The constraint here means:
#   'LIB1_NORM * LIB1_PAR3 * 55663.9^2 / (1.+5.85)' is the fitted "lib.dust.MBB.SED" L_IR, 
#   and 'LIB2_PAR5' is "lib.Radio.mm.FIR.lines.SED" L_IR, 
#   so their ratio is the scaling of L_IR for the "lib.Radio.mm.FIR.lines.SED". 
#   Then, without considering the scaling of L_IR, we still need to scale dL^2 / (1+z)
#   for the flux densities in "lib.Radio.mm.FIR.lines.SED", because at a given L_IR, 
#   if the source is at higher z and larger dL, the flux densities will decrease, 
#   so we need to mulitply 'LIB2_PAR4^2 / 55663.9^2 * (1.+5.85)', 
#   where 'LIB2_PAR4' is the model dL at z=LIB2_PAR3 and '55663.9' is the source dL. 


~/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_plot_SED_fitting_results_for_michi2_v05.py \
"fit_with_CO.out" \
-out "results_fit_with_CO/fit_with_CO.pdf"


~/Cloud/Github/Crab.Toolkit.michi2/bin/michi2_plot_SED_fitting_results_for_michi2_v05.py \
"fit_with_CO.out" \
-out "results_fit_with_CO/fit_with_CO.best.pdf" -only-best


mv best-* "results_fit_with_CO/"




