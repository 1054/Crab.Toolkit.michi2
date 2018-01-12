#!/bin/bash
# 

if [[ $(uname) == "Darwin" ]]; then
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_compute_delta_chisq_1sigma/michi2_compute_delta_chisq_1sigma_mac" bin_mac/
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_read_lib_SED/michi2_read_lib_SED_mac" bin_mac/
fi


mkdir -p ../lib/python/crabtable/ 2>/dev/null
mkdir -p ../lib/python/crabplot/ 2>/dev/null
cp "$HOME/Softwares/Python/lib/crab/crabtable/CrabTable.py" ../lib/python/crabtable/
cp "$HOME/Softwares/Python/lib/crab/crabplot/CrabPlot.py" ../lib/python/crabplot/


wget -q https://raw.githubusercontent.com/1054/DeepFields.SuperDeblending/master/Softwares/ds9_mac/lumdist_mac -O bin_mac/lumdist_mac
wget -q https://raw.githubusercontent.com/1054/DeepFields.SuperDeblending/master/Softwares/ds9_linux_Glibc_2.14/lumdist_linux_x86_64 -O bin_linux_glibc_2_14/lumdist_linux_x86_64
cp michi2_read_lib_SED lumdist
chmod +x lumdist */lumdist*





