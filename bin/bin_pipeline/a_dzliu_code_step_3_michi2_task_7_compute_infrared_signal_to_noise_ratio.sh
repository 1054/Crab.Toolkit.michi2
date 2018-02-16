#!/bin/bash
# 

set -e


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
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))
list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))


# 
# Loop sources
# 
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    
    source_name="${list_of_source_names[i]}"
    
    # check source directory
    if [[ ! -d "$source_name" ]]; then
        echo "Error! \"$source_name\" was not found!"
        exit
    fi
    
    # cd
    echo cd "$source_name/"
    cd "$source_name/"
    
    # run michi2_compute_infrared_signal_to_noise_ratio.py
    michi2_compute_infrared_signal_to_noise_ratio.py extracted_flux.txt SNR.json
    
    # cd back
    cd "../"
    
done





