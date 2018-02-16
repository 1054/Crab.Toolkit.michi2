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
dir_of_photometry="datatable_photometry"
file_of_source_names="$dir_of_photometry/list_of_source_names.txt"
file_of_source_redshifts="$dir_of_photometry/list_of_source_redshifts.txt"
if [[ ! -d "$dir_of_photometry" ]]; then echo "Error! \"$dir_of_photometry\" was not found!"; exit; fi
if [[ ! -f "$file_of_source_names" ]]; then echo "Error! \"$file_of_source_names\" was not found!"; exit; fi
if [[ ! -f "$file_of_source_redshifts" ]]; then echo "Error! \"$file_of_source_redshifts\" was not found!"; exit; fi
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))
list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))


# 
# Make output dir
# 
dir_of_output="SED_fitting_magphys"
if [[ ! -d "$dir_of_output" ]]; then
    mkdir "$dir_of_output"
fi
cp "$file_of_source_names" "$dir_of_output/list_of_source_names.txt"
cp "$file_of_source_redshifts" "$dir_of_output/list_of_source_redshifts.txt"
echo cd "$dir_of_output"
cd "$dir_of_output"


# 
# Loop sources
# 
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    source_name="${list_of_source_names[i]}"
    if [[ ! -d "${source_name}" ]]; then
        mkdir "${source_name}"
    fi
    echo cd "${source_name}"
    cd "${source_name}"
    # 
    #michi2-deploy-files >/dev/null
    #cp ../../$dir_of_photometry/extracted_flux_for_${source_name}.txt extracted_flux.txt
    #cp ../../$dir_of_photometry/extracted_id_ra_dec_zspec_for_${source_name}.txt source_id_ra_dec_zspec.txt
    cp ../../SED_fitting_michi2/${source_name}/extracted_flux.txt .
    michi2_filter_flux_2sigma_no_radio.py extracted_flux.txt flux_obsframe.dat
    # 
    cd "../"
done



cd "../"





