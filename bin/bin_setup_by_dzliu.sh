#!/bin/bash
# 


mkdir bin_mac/ 2>/dev/null
mkdir bin_linux_glibc_2_14/ 2>/dev/null
if [[ $(uname) == "Darwin" ]]; then
    bash -c "cd \"$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_compute_delta_chisq_1sigma/\"; ./do_Compile"
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_compute_delta_chisq_1sigma/michi2_compute_delta_chisq_1sigma_mac" bin_mac/
    
    bash -c "cd \"$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_read_lib_SED/\"; ./do_Compile"
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_read_lib_SED/michi2_read_lib_SED_mac" bin_mac/
    
    bash -c "cd \"$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_v04/\"; ./do_Compile"
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_v04/michi2_v04_mac" bin_mac/
else
    bash -c "cd \"$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_compute_delta_chisq_1sigma/\"; ./do_Compile"
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_compute_delta_chisq_1sigma/michi2_compute_delta_chisq_1sigma_linux_x86_64" bin_linux_glibc_2_14/
    
    bash -c "cd \"$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_read_lib_SED/\"; ./do_Compile"
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_read_lib_SED/michi2_read_lib_SED_linux_x86_64" bin_linux_glibc_2_14/
    
    bash -c "cd \"$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_v04/\"; ./do_Compile"
    cp "$HOME/Cloud/Github/Crab/AstronomyUtility/MiChi2/michi2_v04/michi2_v04_linux_x86_64" bin_linux_glibc_2_14/
fi


mkdir -p ../lib/python/crabtable/ 2>/dev/null
mkdir -p ../lib/python/crabplot/ 2>/dev/null
cp "$HOME/Softwares/Python/lib/crab/crabtable/CrabTable.py" ../lib/python/crabtable/
cp "$HOME/Softwares/Python/lib/crab/crabplot/CrabPlot.py" ../lib/python/crabplot/


wget -q https://raw.githubusercontent.com/1054/DeepFields.SuperDeblending/master/Softwares/ds9_mac/lumdist_mac -O bin_mac/lumdist_mac
wget -q https://raw.githubusercontent.com/1054/DeepFields.SuperDeblending/master/Softwares/ds9_linux_Glibc_2.14/lumdist_linux_x86_64 -O bin_linux_glibc_2_14/lumdist_linux_x86_64

cp bin_portable.sh michi2_v04
cp bin_portable.sh michi2_read_lib_SED
cp bin_portable.sh michi2_compute_delta_chisq_1sigma
cp bin_portable.sh lumdist

chmod +x * */*





