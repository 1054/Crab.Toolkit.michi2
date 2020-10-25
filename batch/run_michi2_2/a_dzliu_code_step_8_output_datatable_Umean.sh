#!/bin/bash
# 
# run this after running
#  ../a_dzliu_code_step_1_get_photometry*.sh
#  ../a_dzliu_code_step_2_extract_photometry*.sh
#  ../a_dzliu_code_step_2_merge_photometry*.sh
#  ../a_dzliu_code_step_3_michi2_task_1_prepare_directories.sh
#  ../a_dzliu_code_step_3_michi2_task_3_run_fitting_on_mac.sh
#  ../a_dzliu_code_step_8_output_datatable_Umean.sh
# 
# 

set -e

data_dir="./SED_fitting_michi2"

output_dir="datatable_Umean"
if [[ ! -d "$output_dir" ]]; then
    mkdir -p "$output_dir"
fi


# 
# set output param list
# 
output_params=(Umin_cold Umean_total fPDR_total LAGN LIR_total LIR_warm LIR_cold Mdust_warm Mdust_cold Mdust_total Mstar)


# 
# read list of source names
# 
list_of_source_names=($(cat "$data_dir/list_of_source_names.txt" | grep -v '^#' | awk '{print $1;}'))
#echo ${list_of_source_names[@]}

for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    
    source_name="${list_of_source_names[i]}"
    new_source_name=$(echo "${list_of_source_names[i]}" | sed -e 's/ID_/V20-ID/g')
    
    #if [[ "${source_name}" == "ID_818" ]] || \
    #   [[ "${source_name}" == "ID_25015" ]] || \
    #   [[ "${source_name}" == "ID_27172" ]] || \
    #   [[ "${source_name}" == "ID_51670" ]]; then
    #    echo "Skipping ${source_name} because of poor IR coverage or fitting."
    #    continue
    #fi
    
    echo "${source_name}"
    
    # loop parameters
    for output_param in ${output_params[@]}; do
        
        data_file="$data_dir/${list_of_source_names[i]}/results_fit_5/best-fit_param_${output_param}.txt"
        
        output_file="${output_dir}/output_datatable_${output_param}.txt"
        
        # check file existence
        if [[ ! -f "$data_file" ]]; then
            echo "Error! \"$data_file\" was not found!"
            exit 255
        fi
        
        # print column header
        if [[ $i == 0 ]]; then
            
            # bakcup previous file
            if [[ -f "${output_file}" ]]; then
                mv "${output_file}" "${output_file}.backup"
            fi
            
            # write output file
            header_str=$(printf "%-15s" "Source")
            cat $data_file | head -n 1 | tail -n 1 | sed -e 's/|/ /g' | sed -e "s/^/# $header_str /g" >  "${output_file}"
            
        fi
        
        # concatenate data to output file
        header_str=$(printf "%-15s" "${new_source_name}")
        cat $data_file | head -n 5 | tail -n 1 | sed -e "s/^/  $header_str /g" >> "${output_file}"
        
    done
    
done






# 
# write readme file
# 
cat << EOF > "${output_dir}/readme.txt"
fitted with the scripts
  ../a_dzliu_code_step_1_get_photometry*.sh
  ../a_dzliu_code_step_2_extract_photometry*.sh
  ../a_dzliu_code_step_2_merge_photometry*.sh
  ../a_dzliu_code_step_3_michi2_task_1_prepare_directories.sh
  ../a_dzliu_code_step_3_michi2_task_3_run_fitting_on_mac.sh
  ../a_dzliu_code_step_8_output_datatable_Umean.sh

date $(date +"%Y%m%d %Hh%Mm%Ss %Z")

EOF



# 
# Done
# 
echo "Output to \"output_datatable_*.txt\"!"
#open "${output_file}"




