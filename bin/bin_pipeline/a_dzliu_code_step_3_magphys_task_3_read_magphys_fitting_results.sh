#!/bin/bash
# 

if [[ ! -d "SED_fitting_magphys_20180110" ]]; then
    echo "Error! \"SED_fitting_magphys_20180110\" was not found!"; exit
fi

cd "SED_fitting_magphys_20180110"


file_of_source_names="list_of_source_names.txt"
file_of_source_redshifts="list_of_source_redshifts.txt"
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))
list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))


for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    source_name=$(basename "${list_of_source_names[i]}")
    if [[ ! -d "${source_name}" ]]; then
        echo "Error! Directory \"${source_name}\" was not found!"; exit 1
    fi
    # 
    #if [[ "${source_name}" != "PACS-164" ]]; then
    #    continue
    #fi
    # 
    echo cd "${source_name}"
    cd "${source_name}/"
    # 
    if [[ ! -f "magphys_fitting/fit_1_with_flux_obsframe/1.fit" ]]; then
        echo "Error! File \"magphys_fitting/fit_1_with_flux_obsframe/1.fit\" was not found!"; exit 1
    fi
    # 
    # Minimum chi2 best fit
    source_name_str=$(printf "%-15s" "$source_name")
    source_header=$(printf "%-13s" "Source")
    cat magphys_fitting/fit_1_with_flux_obsframe/1.sed | grep -A6 "Main parameters of this model" | tail -n 5 | grep -v "^#$" | grep    "^#"  | sed -e 's/\./ /g' | perl -p -e 's/[^a-zA-Z0-9_ ]/_/g' | sed -e 's/^_/#/g' | awk '{print}' ORS=' ' | sed -e 's/$/ /' >  best-fit_param_misc.txt
    cat magphys_fitting/fit_1_with_flux_obsframe/1.sed | grep -A6 "Main parameters of this model" | tail -n 5 | grep -v "^#$" | grep -v "^#"                                                                              | awk '{print}' ORS=' ' | sed -e 's/$/ /' >> best-fit_param_misc.txt
    # 
    # lgMstar P2p5 P16 P50 P84 P97p5
    printf "# %s %s %s %s %s\n" "lgMstar_P2p5" "lgMstar_P16" "lgMstar_P50" "lgMstar_P84" "lgMstar_P97p5"   >  best-fit_param_lgMstar.txt
    cat magphys_fitting/fit_1_with_flux_obsframe/1.fit | head -n 383 | tail -n 1                           >> best-fit_param_lgMstar.txt
    # 
    # lgMdust P2p5 P16 P50 P84 P97p5
    printf "# %s %s %s %s %s\n" "lgMdust_P2p5" "lgMdust_P16" "lgMdust_P50" "lgMdust_P84" "lgMdust_P97p5"   >  best-fit_param_lgMdust.txt
    cat magphys_fitting/fit_1_with_flux_obsframe/1.fit | head -n 790 | tail -n 1                           >> best-fit_param_lgMdust.txt
    # 
    # lgLdust P2p5 P16 P50 P84 P97p5
    printf "# %s %s %s %s %s\n" "lgLdust_P2p5" "lgLdust_P16" "lgLdust_P50" "lgLdust_P84" "lgLdust_P97p5"   >  best-fit_param_lgLdust.txt
    cat magphys_fitting/fit_1_with_flux_obsframe/1.fit | head -n 446 | tail -n 1                           >> best-fit_param_lgLdust.txt
    # 
    # Tdust P2p5 P16 P50 P84 P97p5
    printf "# %s %s %s %s %s\n" "Tdust_P2p5" "Tdust_P16" "Tdust_P50" "Tdust_P84" "Tdust_P97p5"             >  best-fit_param_Tdust.txt
    cat magphys_fitting/fit_1_with_flux_obsframe/1.fit | head -n 1212 | tail -n 1                          >> best-fit_param_Tdust.txt
    # 
    # lgSFR P2p5 P16 P50 P84 P97p5
    printf "# %s %s %s %s %s\n" "lgSFR_P2p5" "lgSFR_P16" "lgSFR_P50" "lgSFR_P84" "lgSFR_P97p5"             >  best-fit_param_lgSFR.txt
    cat magphys_fitting/fit_1_with_flux_obsframe/1.fit | head -n 853 | tail -n 1                           >> best-fit_param_lgSFR.txt
    # 
    # best-fit SED
    pdfcrop -margin 5 magphys_fitting/fit_1_with_flux_obsframe/1.pdf best-fit_SED.pdf
    # 
    # 
    cd "../"
    # 
done







