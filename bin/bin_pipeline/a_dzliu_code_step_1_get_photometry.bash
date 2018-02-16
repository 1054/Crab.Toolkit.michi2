#!/bin/bash
# 

if [[ ! -f "Magellan_SBs_fluxes_for_SED_fitting_mJy.fits" ]]; then
    echo "Error! \"Magellan_SBs_fluxes_for_SED_fitting_mJy.fits\" was not found!"
    exit
fi


if [[ ! -d SED_fitting_michi2 ]]; then mkdir SED_fitting_michi2; fi

cd SED_fitting_michi2


source ~/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash


list_of_source=($(cat "../Magellan_SBs_fluxes_for_SED_fitting.txt" | grep -v "^#" | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 1)) # id_Laigle
#list_of_redshift=($(cat "../Magellan_SBs_fluxes_for_SED_fitting.txt" | grep -v "^#" | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 4)) # 
list_of_redshift=($(cat "../list_of_redshifts.txt" | grep -v "^#" | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 1)) # 


printf "# %-25s\n" "Source" > "list_of_source_names.txt"
printf "# %-15s   %-s\n" "z_spec" "Source" > "list_of_source_redshifts.txt"

for (( i=0; i<${#list_of_source[@]}; i++ )); do
    
    if [[ ! -d "ID_${list_of_source[i]}" ]]; then
        mkdir "ID_${list_of_source[i]}"
    fi
    
    if [[ ! -f "ID_${list_of_source[i]}/extracted_flux.txt" ]]; then
        michi2-extract-flux -cat "../Magellan_SBs_fluxes_for_SED_fitting_mJy.fits" \
                            -id "{'id_Laigle': ${list_of_source[i]}, 'NUMBER': ${list_of_source[i]}}" \
                            -out "ID_${list_of_source[i]}"
        mv "ID_${list_of_source[i]}/extracted_flux.txt" "ID_${list_of_source[i]}/extracted_flux_v1.txt"
        cat "ID_${list_of_source[i]}/extracted_flux_v1.txt" | grep -v "VLA_1.4_GHz_2" > "ID_${list_of_source[i]}/extracted_flux.txt"
    fi
    
    #cp "../../SED_fitting_spdb/fit_matrix_HDFN/fit_sed_data_detected_${list_of_source[i]}.txt" .
    #cp "../../SED_fitting_spdb/fit_matrix_HDFN/fit_sed_data_undetect_${list_of_source[i]}.txt" .
    
    printf "  %-25s\n" "ID_${list_of_source[i]}" >> "list_of_source_names.txt"
    printf "  %-15.4f   %-s\n" "${list_of_redshift[i]}" "ID_${list_of_source[i]}" >> "list_of_source_redshifts.txt"
    
    #break
    
done

