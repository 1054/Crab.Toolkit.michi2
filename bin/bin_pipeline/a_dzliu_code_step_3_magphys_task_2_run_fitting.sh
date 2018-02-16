#!/bin/bash
# 

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
dir_of_output="SED_fitting_magphys"
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
echo cd "$dir_of_output"
cd "$dir_of_output"


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
    # 
    # 
    $HOME/Softwares/magphys/magphys_lowz_go flux_obsframe.dat -redshift ${list_of_source_redshifts[i]}
    #
    # 
    cd "../"
done



cd "../"





