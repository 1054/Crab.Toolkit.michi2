#!/bin/bash
#

./michi2_v03 -obs "flux_co.dat" \
             -lib "lib_z_1.500_co_dvddr_5.0_dVDoW_50.0_Faster.lvg" \
                  "lib_z_1.500_co_dvddr_5.0_dVDoW_50.0_Faster.lvg" \
             -out "fit_two_components.out" \
             -constraint LIB1 PAR1 LT LIB2 PAR1
