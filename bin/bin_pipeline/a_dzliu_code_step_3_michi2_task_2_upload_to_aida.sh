#!/bin/bash
# 

#ssh -t "149.217.40.110" \
#        "mkdir -p /disk1/$USER/Works/AlmaCosmos/SED_Fitting/Sample_SpireLines_R52/SED_fitting_michi2/"


rsync -avz --stats --progress -e "ssh" \
        --include 'list_of_source_names.txt' \
        --include 'list_of_source_redshifts.txt' \
        --exclude '*/*/' \
        --include '*/' \
        --include '*/extracted_flux.txt' \
        --include '*/flux_obsframe.dat' \
        --include '*/source_id_ra_dec_zspec.txt' \
        --include '*/fit_5.out' \
        --exclude '*' \
        "SED_fitting_michi2/" \
        "149.217.40.110:/disk1/$USER/Works/AlmaCosmos/SED_Fitting/Sample_SpireLines_R52/SED_fitting_michi2/" # --dry-run # --verbose


#scp "$HOME/Cloud/Github/Crab.Toolkit.michi2/bin/bin_aida/a_dzliu_code_for_running_SED_fitting_michi2_with_SiebenmorgenAGN_on_aida.sh" \
#        "149.217.40.110:/disk1/$USER/Works/AlmaCosmos/SED_Fitting/Sample_SpireLines_R52/SED_fitting_michi2/"


echo "Then please run \"SED_fitting_michi2/a_dzliu_code_for_running_SED_fitting_michi2_with_SiebenmorgenAGN_on_aida.sh\" on aida40110 machine!"

