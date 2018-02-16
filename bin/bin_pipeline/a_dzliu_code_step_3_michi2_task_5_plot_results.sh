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


# 
# Loop sources
# 
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    source_name="${list_of_source_names[i]}"
    source_redshift="${list_of_source_redshifts[i]}"
    if [[ ! -d "${source_name}" ]]; then
        echo "Error! Directory \"${source_name}\" was not found!"; exit 1
    fi
    # 
    echo cd "${source_name}"
    cd "${source_name}"
    # 
    if [[ ! -f "extracted_flux.txt" ]]; then 
        echo "Error! \"extracted_flux.txt\" was not found under directory \"${source_name}\"!"; exit 1
    fi
    #if [[ ! -f "flux_obsframe.dat" ]]; then 
    #    echo "Error! \"flux_obsframe.dat\" was not found under directory \"${source_name}\"!"; exit 1
    #fi
    if [[ ! -f "fit_5.out" ]]; then 
        echo "Warning! \"fit_5.out\" was not found under directory \"${source_name}\"! Skip and continue!"
        cd "../"
        continue
    fi
    # 
    yrange=()
    if [[ $(awk "BEGIN {if(($source_redshift)<0.1) print 1; else print 0;}") -eq 1 ]]; then
        yrange=("-yrange" "1e-3" "1e7")
    fi
    # 
    if [[ "$*" == *"-overwrite"* ]]; then
        rm -rf best-fit* Plot* obj*
    fi
    # 
    if [[ ! -f "fit_5.pdf" ]]; then
        michi2-plot-fitting-results fit_5.out -flux extracted_flux.txt -source "${source_name}" ${yrange[@]}
    fi
    # 
    if [[ ! -f "fit_5.best.pdf" ]]; then
        michi2-plot-fitting-results fit_5.out -out fit_5.best.pdf -only-best -flux extracted_flux.txt -source "${source_name}" ${yrange[@]}
    fi
    # 
    cd "../"
done



cd "../"





