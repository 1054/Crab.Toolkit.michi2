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
# Check necessary files
# 
dir_of_output="."
if [[ ! -d "$dir_of_output" ]]; then
    echo "Error! Directory \"$dir_of_output\" was not found!"; exit 1
fi
file_of_source_names="$dir_of_output/list_of_source_names.txt"
file_of_source_redshifts="$dir_of_output/list_of_source_redshifts.txt"
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))
list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))


# 
# cd output dir
# 
#echo cd "$dir_of_output"
#cd "$dir_of_output"


# 
# Loop sources
# 
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    source_name="${list_of_source_names[i]}"
    if [[ ! -d "${source_name}" ]]; then
        echo "Error! Directory \"${source_name}\" was not found!"; exit 1
    fi
    # 
    #if [[ "${source_name}" != "PACS-164" ]]; then
    #    continue
    #fi
    # 
    echo cd "${source_name}"
    cd "${source_name}"
    # 
    if [[ ! -f "extracted_flux.txt" ]]; then 
        echo "Error! \"extracted_flux.txt\" was not found under directory \"${source_name}\"!"; exit 1
    fi
    if [[ ! -f "flux_obsframe.dat" ]]; then 
        echo "Error! \"flux_obsframe.dat\" was not found under directory \"${source_name}\"!"; exit 1
    fi
    if [[ ! -f "fit_5.out" ]]; then 
        echo "Error! \"fit_5.out\" was not found under directory \"${source_name}\"!"; cd "../"; continue # exit 1
    fi
    if [[ ! -f "fit_flux_obsframe.dat" ]]; then 
        cp "flux_obsframe.dat" "fit_flux_obsframe.dat"
    fi
    # 
    #michi2_filter_flux_2sigma.py extracted_flux.txt plot_flux_obsframe.dat
    cp extracted_flux.txt plot_flux_obsframe.dat
    cp plot_flux_obsframe.dat flux_obsframe.dat
    # 
    michi2-plot-results-of-fitting-5-components fit_5.out -source "${source_name}"
    # 
    cd "../"
done



cd "../"





