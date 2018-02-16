#!/bin/bash
# 

set -e


# 
# Check software dependencies
# 
if [[ ! -f "$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash" ]]; then
    echo "Error! \"$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash\" was not found! Please clone from \"https://github.com/1054/Crab.Toolkit.michi2\"!"
    exit
fi
source "$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash"


# 
# 
# 
echo "See \"a_dzliu_code_step_2_extract_photometry_data.py\" and \"a_dzliu_code_step_2_task_2_get_SpireLines_photometry_RUN_THIS_MANUALLY.txt\""

echo "Now preparing \"list_of_source_names.txt\" and \"list_of_source_redshifts.txt\""

#find SED_fitting_michi2 -type d -mindepth 1 -maxdepth 1 > SED_fitting_michi2/list_of_source_names.txt
#find SED_fitting_michi2 -type d -mindepth 1 -maxdepth 1 > SED_fitting_michi2/list_of_source_names.txt

echo "# Source"             >  SED_fitting_michi2/list_of_source_names.txt
echo "Arp193"               >> SED_fitting_michi2/list_of_source_names.txt
echo "Arp220"               >> SED_fitting_michi2/list_of_source_names.txt
echo "IRASF17207-0014"      >> SED_fitting_michi2/list_of_source_names.txt
echo "NGC1614"              >> SED_fitting_michi2/list_of_source_names.txt
echo "NGC2623"              >> SED_fitting_michi2/list_of_source_names.txt
echo "NGC6240"              >> SED_fitting_michi2/list_of_source_names.txt
echo "UGC05101"             >> SED_fitting_michi2/list_of_source_names.txt

echo "# Source"             >  SED_fitting_michi2/list_of_source_redshifts.txt
echo "0.023299"             >> SED_fitting_michi2/list_of_source_redshifts.txt
echo "0.018126"             >> SED_fitting_michi2/list_of_source_redshifts.txt
echo "0.042810"             >> SED_fitting_michi2/list_of_source_redshifts.txt
echo "0.015938"             >> SED_fitting_michi2/list_of_source_redshifts.txt
echo "0.018509"             >> SED_fitting_michi2/list_of_source_redshifts.txt
echo "0.024480"             >> SED_fitting_michi2/list_of_source_redshifts.txt
echo "0.039367"             >> SED_fitting_michi2/list_of_source_redshifts.txt



