#!/bin/bash
# 


rsync -ravz --prune-empty-dirs --stats --progress -e "ssh" \
        --exclude 'list_of_source_names.txt' \
        --exclude 'list_of_source_redshifts.txt' \
        --exclude 'a_dzliu_code_*' \
        --exclude 'BACKUP*' \
        --exclude '.*' \
        --include '*' \
        --exclude '*' \
        "149.217.40.110:/disk1/$USER/Works/AlmaCosmos/SED_Fitting/Sample_SpireLines_R52/SED_fitting_michi2/" \
        "SED_fitting_michi2/" --dry-run


