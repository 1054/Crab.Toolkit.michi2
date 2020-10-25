#!/bin/bash
# 

if [[ ! -d SED_fitting_michi2 ]]; then
    echo "Error! Directory not found: \"SED_fitting_michi2\""
    exit 255
fi
echo cd SED_fitting_michi2
cd SED_fitting_michi2

set -e


# 
# Check software dependencies
# 
Github_dir="$HOME/Cloud/Github"
if [[ $(hostname) == isaac* ]]; then
    Github_dir="/u/$USER/Cloud/Github"
fi
# 
if [[ ! -f "$Github_dir/Crab.Toolkit.michi2/SETUP.bash" ]]; then
    echo "Error! \"$Github_dir/Crab.Toolkit.michi2/SETUP.bash\" was not found! Please clone from \"https://github.com/1054/Crab.Toolkit.michi2\"!"
    exit
else
    source "$Github_dir/Crab.Toolkit.michi2/SETUP.bash"
fi


# 
# Check necessary files
# 
file_of_source_names="list_of_source_names.txt"
file_of_source_redshifts="list_of_source_redshifts.txt"
if [[ ! -f "$file_of_source_names" ]]; then echo "Error! \"$file_of_source_names\" was not found!"; exit; fi
if [[ ! -f "$file_of_source_redshifts" ]]; then echo "Error! \"$file_of_source_redshifts\" was not found!"; exit; fi
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | tr -s ' ' | cut -d ' '  -f 1))
list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | tr -s ' ' | cut -d ' '  -f 1))


# 
# Loop sources
# 
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    
    itask=$((i+1))
    
    if [[ x"${list_of_source_names[i]}" == x"" ]]; then
        continue
    fi
    
    #if [[ x"${list_of_source_names[i]}" == x"ID_818" ]]; then
    #    continue # no IR, only Laigle+2016 data
    #fi
    
    # print progress
    echo "Processing \"${list_of_source_names[i]}\"  (${itask} / ${#list_of_source_names[@]})"
    
    # check slurm task ID, only run the number which equals current slurm task ID.
    if [[ x"$SLURM_ARRAY_TASK_ID" != x"" ]]; then
        if [[ $itask -ne $SLURM_ARRAY_TASK_ID ]]; then
            continue
        fi
    fi
    
    # check source dir
    if [[ ! -d "${list_of_source_names[i]}" ]]; then
        echo "Warning! \"${list_of_source_names[i]}\" was not found! Will skip this source!"
        continue
    fi
    
    # cd source dir
    cd "${list_of_source_names[i]}"
    
    
    # prepare flux_obsframe.dat
    if [[ ! -f "flux_obsframe.dat" ]]; then
        michi2_filter_flux_2sigma_fit_infrared_upper_limits.py datatable_photometry.txt flux_obsframe.dat
    fi
    
    
    # run SED fitting
    if [[ ! -f "fit_5.out" ]]; then
        
        #michi2-deploy-files
        
        #michi2_filter_flux_2sigma_fit_infrared_upper_limits.py datatable_photometry.txt fit_5.in
        
        michi2-run-SED-fitting-v5 \
            -redshift ${list_of_source_redshifts[i]} \
            -flux flux_obsframe.dat \
            -parallel 2 \
            -obj-name "${list_of_source_names[i]}" \
            -lib-stellar BC03.MultiAge \
            -sampling 15000
            
        
        sleep 2
          
          
        # delete this figure so as to remake figure with "datatable_photometry.txt"
        if [[ -f "results_fit_5/fit_5.pdf" ]]; then
            rm "results_fit_5/fit_5.pdf"
        fi
        if [[ -f "results_fit_5/fit_5.best.pdf" ]]; then
            rm "results_fit_5/fit_5.best.pdf"
        fi
        
    else
        
        echo "Found existing \"fit_5.out\" under directory \"${list_of_source_names[i]}\"! Will skip this source!"
        
    fi
    
    # remake figures
    #rm "results_fit_5/fit_5.pdf"
    if [[ ! -f "results_fit_5/fit_5.pdf" ]]; then
        michi2_plot_SED_fitting_results_for_michi2_v05.py fit_5.out -flux datatable_photometry.txt -source "${list_of_source_names[i]}" -out results_fit_5/fit_5.pdf
        sleep 1
    else
        echo "Found existing \"fit_5.pdf\" under directory \"${list_of_source_names[i]}\"!"
    fi
    
    #rm "results_fit_5/fit_5.best.pdf"
    if [[ ! -f "results_fit_5/fit_5.best.pdf" ]]; then
        michi2_plot_SED_fitting_results_for_michi2_v05.py fit_5.out -flux datatable_photometry.txt -source "${list_of_source_names[i]}" -out results_fit_5/fit_5.best.pdf -only-best
        sleep 1
    else
        echo "Found existing \"fit_5.best.pdf\" under directory \"${list_of_source_names[i]}\"!"
    fi
    
    # calc SNR_FIRMM
    #if [[ -f "SNR_FIRMM.txt" ]]; then
    #    rm "SNR_FIRMM.txt"
    #fi
    if [[ ! -f "SNR_FIRMM.txt" ]]; then
        michi2_compute_infrared_signal_to_noise_ratio.py datatable_photometry.txt SNR_FIRMM.txt
    else
        echo "Found existing \"SNR_FIRMM.txt\" under directory \"${list_of_source_names[i]}\"!"
    fi
    
    # cd back
    cd "../"
    
    #break
    
done