#!/bin/bash
# 

source ~/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash

cd ID_628146/

mkdir SED_fitting_michi2/
cd SED_fitting_michi2/

cp ../datatable_id_ra_dec_zspec.txt ../datatable_photometry.txt .

redshift=$(grep -v '^#' "datatable_id_ra_dec_zspec.txt" | awk '{ print $4 }')

michi2-run-SED-fitting-v5 -flux datatable_photometry.txt -redshift $redshift -trial

michi2-run-SED-fitting-v5 -flux datatable_photometry.txt -redshift $redshift -parallel 2 -sampling 3000 -stellar BC03.400Myr -overwrite


