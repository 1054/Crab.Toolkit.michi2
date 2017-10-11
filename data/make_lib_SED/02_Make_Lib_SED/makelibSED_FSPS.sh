#!/bin/bash
# 

tar -xzf ../01_SED_Models/FSPS/fsps_ssp/SSP_Padova_BaSeL_Chabrier_allZ.tar.gz -C FSPSspec

idl -e 'makelibSED_FSPS'

