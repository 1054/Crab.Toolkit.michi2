#!/bin/bash
# 


rsync -ravz --prune-empty-dirs --stats --progress -e "ssh" \
        --include 'SED_fitting_michi2/' \
        --include 'SED_fitting_michi2/ID_*/' \
        --include 'SED_fitting_michi2/ID_*/fit_*' \
        --include 'SED_fitting_michi2/ID_*/Plot_*.pdf' \
        --include 'SED_fitting_michi2/ID_*/best-fit_*' \
        --include 'SED_fitting_michi2/ID_*/obj_*/' \
        --include 'SED_fitting_michi2/ID_*/obj_*/*' \
        --exclude '*' \
        "149.217.40.110:/disk1/$USER/Works/DeepFields/Works_co/2017_Emanuele_and_Antonello/20180215_SED_fitting_michi2_evolving_qIR/" \
        "./" --dry-run


