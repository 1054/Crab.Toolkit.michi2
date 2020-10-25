#!/bin/bash
# 

if [[ ! -d SED_fitting_michi2 ]]; then
    mkdir SED_fitting_michi2
fi
echo cd SED_fitting_michi2
cd SED_fitting_michi2

set -e

overwrite=0

input_cat="datatable_photometry/datatable_crossmatched_final.fits"

if [[ ! -f "extracted_flux.log" ]] || [[ $overwrite -ge 1 ]]; then
     ~/Cloud/Github/Crab.Toolkit.michi2/bin/michi2-extract-flux -catalog "../$input_cat"
fi
if [[ ! -f "extracted_wavelength.txt" ]] || [[ $overwrite -ge 1 ]]; then
     topcat -stilts tpipe in="../$input_cat" cmd="keepcols \"WAVELENGTH_ALMA\"" out="extracted_wavelength.txt" ofmt=ascii
fi
if [[ ! -f "extracted_zspec.txt" ]] || [[ $overwrite -ge 1 ]]; then
     topcat -stilts tpipe in="../$input_cat" cmd="keepcols \"z ID\"" out="extracted_zspec.txt" ofmt=CSV
fi
if [[ ! -f "extracted_ra.txt" ]] || [[ $overwrite -ge 1 ]]; then
     topcat -stilts tpipe in="../$input_cat" cmd="keepcols \"RA\"" out="extracted_ra.txt" ofmt=ascii
fi
if [[ ! -f "extracted_dec.txt" ]] || [[ $overwrite -ge 1 ]]; then
     topcat -stilts tpipe in="../$input_cat" cmd="keepcols \"Dec\"" out="extracted_dec.txt" ofmt=ascii
fi



if [[ -f "list_of_source_names.txt" ]] && [[ $overwrite -ge 1 ]]; then
     mv "list_of_source_names.txt" "list_of_source_names.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
     mv "list_of_source_redshifts.txt" "list_of_source_redshifts.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
     mv "list_of_source_zspec.txt" "list_of_source_zspec.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
     mv "list_of_source_radec.txt" "list_of_source_radec.txt.backup."$(date "+%Y%m%d.%H%M%S.%Z")
fi


if [[ ! -f "list_of_source_names.txt" ]]; then
     echo "# SourceName" > "list_of_source_names.txt"
     echo "# z   SourceName" > "list_of_source_redshifts.txt"
     echo "# zspec   SourceName" > "list_of_source_zspec.txt"
     echo "# RA   Dec   SourceName" > "list_of_source_radec.txt"
     i=1
     while [[ -f "extracted_flux_for_obj_at_row_${i}.txt" ]]; do
          ID_Master=$(cat "extracted_flux_for_obj_at_row_${i}.info" | grep "^ID = " | sed -e 's/ID = //g')
          if [[ $(awk "BEGIN {if ($ID_Master>0) print 1; else print 0;}") -eq 1 ]]; then
               ALMA_Wavelength=$(printf "%-20.6f" $(cat "extracted_wavelength.txt" | grep -v "^#" | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 1))
               Source_ra=$(cat "extracted_ra.txt" | grep -v "^#" | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 1)
               Source_dec=$(cat "extracted_dec.txt" | grep -v "^#" | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 1)
               Source_zspec=$(cat "extracted_zspec.txt" | tail -n +2 | head -n $i | tail -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ',' -f 1); if [[ x$(echo "${Source_zspec}" | sed -e 's/[^0-9.+-eE]//g') == x ]]; then Source_zspec="-99"; fi; if [[ $(awk "BEGIN {if(${Source_zspec}>=10.0) print 1; else print 0;}") -eq 1 ]]; then Source_zspec="-100"; fi
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
                    echo "${Source_zspec}   ID_${ID_Master}" >> "list_of_source_redshifts.txt"
                    echo "${Source_zspec}   ID_${ID_Master}" >> "list_of_source_zspec.txt"
                    echo "${Source_ra}   ${Source_dec}   ID_${ID_Master}" >> "list_of_source_radec.txt"
                    echo "# ID  RA  Dec  zspec" > "ID_${ID_Master}/datatable_id_ra_dec_zspec.txt"
                    echo "${ID_Master}  ${Source_ra}  ${Source_dec}  ${Source_zspec}" >> "ID_${ID_Master}/datatable_id_ra_dec_zspec.txt"
               else
                    dupl=$((dupl+1))
                    echo "ID_${ID_Master} (dupl $dupl)"
               fi
               if [[ ! -f "ID_${ID_Master}/datatable_photometry.txt" ]]; then
                    head -n 1 "extracted_flux_for_obj_at_row_${i}.txt" > "ID_${ID_Master}/datatable_photometry.txt"
               fi
               # set A3COSMOS photometry data points
               cat "extracted_flux_for_obj_at_row_${i}.txt" | grep -v "^#" | perl -p -e "s/([ ]+)nan([ ]+)([0-9.+-eE]+)([ ]+)([0-9.+-eE]+)([ ]+)(mJy)([ ]+)unknown_input_str_FLUX_ALMA/  $ALMA_Wavelength \3\4\5\6\7\8Input_cat_row_${i}/g" >> "ID_${ID_Master}/datatable_photometry.txt"
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

echo "Done!"


