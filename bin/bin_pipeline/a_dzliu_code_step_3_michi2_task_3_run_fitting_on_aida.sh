#!/bin/bash
# 


# 
# Define output directory
# 
dir_of_output="SED_fitting_michi2"

if [[ ! -d "$dir_of_output" ]]; then
    echo "Error! Directory \"$dir_of_output\" was not found!"; exit 1
fi

echo cd "$dir_of_output"
cd "$dir_of_output"




# 
# Check software dependencies
# 
if [[ ! -f "$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash" ]]; then
    echo "Error! \"$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash\" was not found! Please clone from \"https://github.com/1054/Crab.Toolkit.michi2\"!"
    exit
fi
source "$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash"


# 
# Read source names and redshifts
# 
file_of_source_names="list_of_source_names.txt"
if [[ ! -f "$file_of_source_names" ]]; then
    echo "Error! File \"$file_of_source_names\" was not found under $(pwd)!"; exit 1
fi
file_of_source_redshifts="list_of_source_redshifts.txt"
if [[ ! -f "$file_of_source_redshifts" ]]; then
    echo "Error! File \"$file_of_source_redshifts\" was not found under $(pwd)!"; exit 1
fi
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | tr -s ' ' | cut -d ' ' -f 1))
list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | tr -s ' ' | cut -d ' ' -f 1))



for (( i=0; i<${#list_of_source_names[@]}; i++ )); do
    
    if [[ ! -f "fit_5.out" ]]; then
        # 
        echo 
        echo 
        echo "***********************************************"
        echo "cd \"ID_${list_of_source_names[i]}/\""
        echo "***********************************************"
        cd "ID_${list_of_source_names[i]}/"
        # 
        echo "michi2_filter_flux_2sigma.py extracted_flux.txt fit_5.in"
        michi2_filter_flux_2sigma.py extracted_flux.txt fit_5.in
        # 
        echo "michi2-run-fitting-5-components -redshift ${list_of_source_redshifts[i]} -parallel 12"
        michi2-run-fitting-5-components -redshift ${list_of_source_redshifts[i]} -parallel 12
        # 
        echo "michi2-plot-fitting-results fit_5.out -source \"${list_of_source_names[i]}\""
        michi2-plot-fitting-results fit_5.out -source "${list_of_source_names[i]}"
        # 
        echo "cd \"../\""
        cd "../"
    fi
    
    #break
    
done


