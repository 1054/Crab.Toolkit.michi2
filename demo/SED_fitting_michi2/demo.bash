#!/bin/bash
#

source ~/Cloud/Github/Crab.toolkit.michi2/SETUP.bash

cd ID_628146/

mkdir -p SED_fitting_michi2/
cd SED_fitting_michi2/

cp ../datatable_id_ra_dec_zspec.txt ../datatable_photometry.txt .

redshift=$(grep -v '^#' "datatable_id_ra_dec_zspec.txt" | awk '{ print $4 }')

# Deploy the model libraries (required in a fresh directory — the wrapper's
# auto-deploy does not trigger when none of the libraries are present yet).
michi2-deploy-files SED

michi2-run-SED-fitting-v5 -flux datatable_photometry.txt -redshift $redshift -trial

michi2-run-SED-fitting-v5 -flux datatable_photometry.txt -redshift $redshift -parallel 2 -sampling 3000 -lib-stellar BC03.MultiAge -overwrite
