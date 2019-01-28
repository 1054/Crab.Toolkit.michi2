#!/bin/bash

# check files
if [[ ! -f 'datatable_photometry_with_NOEMA_1mm_with_optical.fits' ]]; then
    echo "Please run step_1 first!"
    exit
fi


source ~/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash


input_cat="datatable_photometry_with_NOEMA_1mm_with_optical.fits"
col_id="id_Liu2018"
col_z="z_Liu2018"
col_ref_z="ref_z_Liu2018"
col_ra="ra_Liu2018"
col_dec="de_Liu2018"


overwrite=1


if [[ ! -d "Multi-wavelength_SEDs" ]]; then 
    mkdir "Multi-wavelength_SEDs"
fi

cd "Multi-wavelength_SEDs"

if [[ ! -f "extracted_flux.log" ]] || [[ $overwrite -ge 1 ]]; then
     michi2-extract-flux-v2 -catalog "../$input_cat"
fi
if [[ ! -f "extracted_zprior.txt" ]] || [[ $overwrite -ge 1 ]]; then
     topcat -stilts tpipe in="../$input_cat" cmd="keepcols \"$col_z $col_ref_z $col_id\"" out="extracted_zprior.txt" ofmt=CSV
fi
if [[ ! -f "extracted_ra.txt" ]] || [[ $overwrite -ge 1 ]]; then
     topcat -stilts tpipe in="../$input_cat" cmd="keepcols \"$col_ra\"" out="extracted_ra.txt" ofmt=ascii
fi
if [[ ! -f "extracted_dec.txt" ]] || [[ $overwrite -ge 1 ]]; then
     topcat -stilts tpipe in="../$input_cat" cmd="keepcols \"$col_dec\"" out="extracted_dec.txt" ofmt=ascii
fi


if [[ -f "list_of_source_names.txt" ]] && [[ $overwrite -ge 1 ]]; then
     mv "list_of_source_names.txt" "list_of_source_names.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
     mv "list_of_source_zphot.txt" "list_of_source_zphot.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
     mv "list_of_source_zspec.txt" "list_of_source_zspec.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
     mv "list_of_source_radec.txt" "list_of_source_radec.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
fi


if [[ ! -f "list_of_source_names.txt" ]]; then
     echo "# SourceName" > "list_of_source_names.txt"
     echo "# zphot   SourceName" > "list_of_source_zphot.txt"
     echo "# zspec   SourceName" > "list_of_source_zspec.txt"
     echo "# RA   Dec   SourceName" > "list_of_source_radec.txt"
     i=1
     while [[ -f "extracted_flux_for_obj_at_row_${i}.txt" ]]; do
          ID_Master=$(cat "extracted_flux_for_obj_at_row_${i}.info" | grep "^$col_id = " | sed -e "s/$col_id = //g")
          if [[ $(awk "BEGIN {if ($ID_Master>0) print 1; else print 0;}") -eq 1 ]]; then
               Source_ra=$(cat "extracted_ra.txt" | grep -v "^#" | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 1)
               Source_dec=$(cat "extracted_dec.txt" | grep -v "^#" | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 1)
               # 
               #Source_zphot=$(cat "extracted_zphot.txt" | tail -n +2 | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ',' -f 1); if [[ x$(echo "${Source_zphot}" | sed -e 's/[^0-9.+-eE]//g') == x ]]; then Source_zphot="-99"; fi
               #Source_zspec=$(cat "extracted_zspec.txt" | tail -n +2 | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ',' -f 1); if [[ x$(echo "${Source_zspec}" | sed -e 's/[^0-9.+-eE]//g') == x ]]; then Source_zspec="-99"; fi; if [[ $(awk "BEGIN {if(${Source_zspec}>=10.0) print 1; else print 0;}") -eq 1 ]]; then Source_zspec="-100"; fi
               Source_zpriors=($(cat "extracted_zprior.txt" | tail -n +2 | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ',' -f 1))
               Source_refzpriors=($(cat "extracted_zprior.txt" | tail -n +2 | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ',' -f 2))
               # tail -n +2 means skip the first line, which is the header line of CSV format 
               # 
               # check folder
               if [[ ! -d "ID_${ID_Master}" ]]; then 
                    mkdir "ID_${ID_Master}"
               fi
               # 
               # check overwritting
               if [[ -f "ID_${ID_Master}/extracted_flux_for_obj_at_row_${i}.txt" ]]; then
                    if [[ $overwrite -ge 1 ]]; then
                         echo "Overwritting! rm \"ID_${ID_Master}/extracted_flux_for_obj_at_row_\"*\".txt\""
                         rm "ID_${ID_Master}/extracted_flux_for_obj_at_row_"*".txt"
                    else
                         echo "Error! Found existing \"ID_${ID_Master}/extracted_flux_for_obj_at_row_${i}.txt\"! We will not overwrite unless set in the script!"
                         echo "Exit!"
                         exit
                    fi
               fi
               # 
               # check duplicated photometry
               # 
               dupl=$(find "ID_${ID_Master}" -maxdepth 1 -name "extracted_flux_for_obj_at_row_*.txt" | wc -l)
               # 
               if [[ $dupl -eq 0 ]]; then 
                    echo "ID_${ID_Master}"
                    echo "ID_${ID_Master}" >> "list_of_source_names.txt"
                    #echo "${Source_zphot}   ID_${ID_Master}" >> "list_of_source_zphot.txt"
                    #echo "${Source_zspec}   ID_${ID_Master}" >> "list_of_source_zspec.txt"
                    echo "${Source_ra}   ${Source_dec}   ID_${ID_Master}" >> "list_of_source_radec.txt"
                    # z_phot, z_spec (obsolete since 20180723)
                    echo "# ID  RA  Dec  zphot" > "ID_${ID_Master}/datatable_id_ra_dec_zphot.txt"
                    echo "# ID  RA  Dec  zspec" > "ID_${ID_Master}/datatable_id_ra_dec_zspec.txt"
                    echo "${ID_Master}  ${Source_ra}  ${Source_dec}  ${Source_zphot}" >> "ID_${ID_Master}/datatable_id_ra_dec_zphot.txt"
                    echo "${ID_Master}  ${Source_ra}  ${Source_dec}  ${Source_zspec}" >> "ID_${ID_Master}/datatable_id_ra_dec_zspec.txt"
                    # z_prior (20180723)
                    if [[ ${#Source_zpriors[@]} -gt 0 ]]; then
                         printf "# %-13s %15s %15s  %15s   %s\n" "ID" "RA" "Dec" "z_prior" "ref_z_prior" > "ID_${ID_Master}/datatable_id_ra_dec_zprior.txt"
                         for (( j = 0; j < ${#Source_zpriors[@]}; j++ )); do
                              printf "%-15d %15.8f %15.8f  %15g   %s\n" "${ID_Master}" "${Source_ra}" "${Source_dec}" "${Source_zpriors[j]}" "${Source_refzpriors[j]}" >> "ID_${ID_Master}/datatable_id_ra_dec_zprior.txt"
                         done
                    fi
               else
                    dupl=$((dupl+1))
                    echo "ID_${ID_Master} (dupl $dupl)"
               fi
               if [[ ! -f "ID_${ID_Master}/datatable_photometry.txt" ]]; then
                    cp "extracted_flux_for_obj_at_row_${i}.txt" "ID_${ID_Master}/datatable_photometry.txt"
               else
                    cat "extracted_flux_for_obj_at_row_${i}.txt" | grep -v "^#" | perl -p -e "s/(.*)(\w)\s*$/\1\2_row_${i}/g" >> "ID_${ID_Master}/datatable_photometry.txt"
               fi
               mv "extracted_flux_for_obj_at_row_${i}.info" "ID_${ID_Master}/"
               mv "extracted_flux_for_obj_at_row_${i}.txt" "ID_${ID_Master}/"
          else
               
               if [[ ! -d "ID_Not_In_Master_Catalog" ]]; then 
                    mkdir "ID_Not_In_Master_Catalog"
               fi
               mv "extracted_flux_for_obj_at_row_${i}.info" "ID_Not_In_Master_Catalog/"
               mv "extracted_flux_for_obj_at_row_${i}.txt" "ID_Not_In_Master_Catalog/"
               
          fi
          i=$((i+1))
     done
     # 
     #mv "list_of_source_names.txt" "list_of_source_names_unsorted.txt"
     #cat "list_of_source_names_unsorted.txt" | grep -v "^#" | sort -V >  "list_of_source_names.txt"
fi



# Note that datatable_photometry.txt still contains duplicated wavelengths, 
# see "a_dzliu_code_step_3_run_sed_fitting_spdb.bash" for solving this.







