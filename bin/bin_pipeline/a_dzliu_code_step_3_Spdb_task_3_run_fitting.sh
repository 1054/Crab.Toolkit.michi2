#!/bin/bash
# 


# 
# Define output directory
# 
dir_of_output="SED_fitting_Spdb"

if [[ ! -d "$dir_of_output" ]]; then
    echo "Error! Directory \"$dir_of_output\" was not found!"; exit 1
fi

echo cd "$dir_of_output"
cd "$dir_of_output"




# 
# Check software dependencies
# 
if [[ ! -f "$HOME/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP.bash" ]]; then
    echo "Error! \"$HOME/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP.bash\" was not found! Please clone from \"https://github.com/1054/DeepFields.SuperDeblending\"!"
    exit
fi
source "$HOME/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP.bash"


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
    
    if [[ ! -f "${list_of_source_names[i]}/fit_plots_HDFN/Plot_SED_"$(echo ${list_of_source_names[i]} | sed -e 's/ID_//g')".pdf" ]]; then
        # 
        echo 
        echo 
        echo "***********************************************"
        echo "cd \"${list_of_source_names[i]}/\""
        echo "***********************************************"
        mkdir "${list_of_source_names[i]}"
        cd "${list_of_source_names[i]}/"
        # 
        cp "../../SED_fitting_michi2/${list_of_source_names[i]}/extracted_flux.txt" datatable_photometry.txt
        echo "# ID RA Dec z_spec" > datatable_id_ra_dec_zspec.txt
        echo $(echo ${list_of_source_names[i]} | sed -e 's/ID_//g') 0.0 0.0 ${list_of_source_redshifts[i]} >> datatable_id_ra_dec_zspec.txt
        # 
        ~/Cloud/Github/DeepFields.SuperDeblending/Softwares/Galsed_Template/deploy_files
        # 
        echo "macro read convert_datatable_for_sed_fitting.sm convert_datatable_for_sed_fitting" | sm
        # 
        ./do_Galsed
        # 
        echo "cd \"../\""
        cd "../"
    fi
    
    #break
    
done


