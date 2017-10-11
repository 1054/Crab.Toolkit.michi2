#!/bin/bash
# 

cp  '/Users/dzliu/Softwares/BC03/bc03/out/out_constant_SFH_wave_flux_at_200Myr.txt' \
    'datatable_wave_AA_flux_Lsun_per_AA.txt'

idl -e 'makelibSED_BC03'

