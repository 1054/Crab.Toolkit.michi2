#!/bin/bash
# 

#ssh -t "149.217.40.110" \
#        "mkdir -p /disk1/$USER/Works/DeepFields/Works_co/2017_Emanuele_and_Antonello/20180215_SED_fitting_michi2_evolving_qIR/"


rsync -avz --stats --progress -e "ssh" \
        --exclude '.*' \
        --include 'a_dzliu_code_step_3_task_3_run_michi2_on_aida40110.bash' \
        --exclude '*/*/*' \
        --include 'SED_fitting_michi2/' \
        --include 'SED_fitting_michi2/*.txt' \
        --exclude 'SED_fitting_michi2/Results' \
        --include 'SED_fitting_michi2/*/' \
        --include 'SED_fitting_michi2/*/extracted_flux*' \
        --include 'SED_fitting_michi2/*/source_id_ra_dec_z*' \
        --include 'SED_fitting_michi2/*/fit_5.in' \
        --include 'SED_fitting_michi2/*/fit_5.out' \
        --exclude 'SED_fitting_michi2/*' \
        "./" \
        "149.217.40.110:/disk1/$USER/Works/DeepFields/Works_co/2017_Emanuele_and_Antonello/20180215_SED_fitting_michi2_evolving_qIR/" # --dry-run # --verbose


#scp "$HOME/Cloud/Github/Crab.Toolkit.michi2/bin/bin_aida/a_dzliu_code_for_running_SED_fitting_michi2_on_aida.sh" \
#        "149.217.40.110:/disk1/$USER/Works/DeepFields/Works_co/2017_Emanuele_and_Antonello/20180215_SED_fitting_michi2_evolving_qIR/"


#echo "Then please run \"a_dzliu_code_for_running_SED_fitting_michi2_on_aida.sh\" on aida40110 machine!"
echo "Then please run \"a_dzliu_code_step_3_task_3_run_michi2_on_aida40110.bash\" on aida40110 machine!"




