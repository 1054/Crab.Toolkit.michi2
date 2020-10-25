#!/bin/bash
# 
# run this after running
#  a_dzliu_code_step_1_get_photometry*.sh
#  a_dzliu_code_step_2_extract_photometry*.sh
#  a_dzliu_code_step_2_merge_photometry*.sh
#  a_dzliu_code_step_3_michi2_task_1_prepare_directories.sh
#  a_dzliu_code_step_3_michi2_task_3_run_fitting_on_mac.sh
#  a_dzliu_code_step_4_michi2_collect_results.sh
# 
# 

set -e

data_dir="./SED_fitting_michi2"

output_dir="datatable_SED_fitting_michi2"
if [[ ! -d "$output_dir" ]]; then
    mkdir -p "$output_dir"
fi


# 
# set output param list
# 
output_params=(Umin_cold Umean_total fPDR_total LAGN LIR_total LIR_warm LIR_cold Umin_cold Mdust_warm Mdust_cold Mdust_total Mstar)


# 
# read list of source names
# 
list_of_source_names=($(cat "$data_dir/list_of_source_names.txt" | grep -v '^#' | awk '{print $1;}'))
list_of_libraries=()
#echo ${list_of_source_names[@]}

for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    
    echo "${list_of_source_names[i]}"
    
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
        header_str=$(printf "%-15s" "${list_of_source_names[i]}")
        cat $data_file | head -n 5 | tail -n 1 | sed -e "s/^/  $header_str /g" >> "${output_file}"
        
    done
    
    
    if [[ $i == 0 ]]; then
        # record list_of_libraries
        IFS='' list_of_libraries+=($(cat "$data_dir/${list_of_source_names[i]}/fit_5.out.info" | grep 'LIB'))
    fi
    
done



# 
# write readme file
# 
echo "Current date time: $(date +"%Y%m%d %Hh%Mm%Ss %Z")" > "${output_dir}/readme.txt"
echo "Fitted with michi2_v05 with libraries:" >> "${output_dir}/readme.txt"
for (( i = 0; i < ${#list_of_libraries[@]}; i++ )); do
    echo "${list_of_libraries[i]}" >> "${output_dir}/readme.txt"
done



# 
# Done
# 
echo "Output to \"${output_dir}/output_datatable_*.txt\"!"
#open "${output_file}"



# 
# cd ${output_dir}
# 
echo cd "${output_dir}"
cd "${output_dir}"



# 
# make All in one datatable
# 
if [[ -f output_datatable_all_in_one.fits ]]; then
    mv output_datatable_all_in_one.fits output_datatable_all_in_one.fits.backup
fi
topcat -stilts tmatchn nin=11 \
        in1=output_datatable_Mstar.txt ifmt1=ascii values1='Source' suffix1='' \
        in2=output_datatable_LAGN.txt ifmt2=ascii values2='Source' suffix2='_2' \
        in3=output_datatable_LIR_total.txt ifmt3=ascii values3='Source' suffix3='_3' \
        in4=output_datatable_LIR_warm.txt ifmt4=ascii values4='Source' suffix4='_4' \
        in5=output_datatable_LIR_cold.txt ifmt5=ascii values5='Source' suffix5='_5' \
        in6=output_datatable_Mdust_total.txt ifmt6=ascii values6='Source' suffix6='_6' \
        in7=output_datatable_Mdust_warm.txt ifmt7=ascii values7='Source' suffix7='_7' \
        in8=output_datatable_Mdust_cold.txt ifmt8=ascii values8='Source' suffix8='_8' \
        in9=output_datatable_Umean_total.txt ifmt9=ascii values9='Source' suffix9='_9' \
        in10=output_datatable_Umin_cold.txt ifmt10=ascii values10='Source' suffix10='_10' \
        in11=output_datatable_fPDR_total.txt ifmt11=ascii values11='Source' suffix11='_11' \
        \
        matcher=exact iref=1 fixcols=dups \
        \
        icmd1='replacecol "param_median" -name "Mstar_median" "toDouble(param_median)"' \
        icmd2='replacecol "param_median" -name "LAGN_median" "toDouble(param_median)"' \
        icmd3='replacecol "param_median" -name "LIR_median" "toDouble(param_median)"' \
        icmd4='replacecol "param_median" -name "LIR_warm_median" "toDouble(param_median)"' \
        icmd5='replacecol "param_median" -name "LIR_cold_median" "toDouble(param_median)"' \
        icmd6='replacecol "param_median" -name "Mdust_median" "toDouble(param_median)"' \
        icmd7='replacecol "param_median" -name "Mdust_warm_median" "toDouble(param_median)"' \
        icmd8='replacecol "param_median" -name "Mdust_cold_median" "toDouble(param_median)"' \
        icmd9='replacecol "param_median" -name "Umean_median" "toDouble(param_median)"' \
        icmd10='replacecol "param_median" -name "Umin_median" "toDouble(param_median)"' \
        icmd11='replacecol "param_median" -name "fPDR_median" "toDouble(param_median)"' \
        \
        icmd1='replacecol "param_best" -name "Mstar_best" "toDouble(param_best)"' \
        icmd2='replacecol "param_best" -name "LAGN_best" "toDouble(param_best)"' \
        icmd3='replacecol "param_best" -name "LIR_best" "toDouble(param_best)"' \
        icmd4='replacecol "param_best" -name "LIR_warm_best" "toDouble(param_best)"' \
        icmd5='replacecol "param_best" -name "LIR_cold_best" "toDouble(param_best)"' \
        icmd6='replacecol "param_best" -name "Mdust_best" "toDouble(param_best)"' \
        icmd7='replacecol "param_best" -name "Mdust_warm_best" "toDouble(param_best)"' \
        icmd8='replacecol "param_best" -name "Mdust_cold_best" "toDouble(param_best)"' \
        icmd9='replacecol "param_best" -name "Umean_best" "toDouble(param_best)"' \
        icmd10='replacecol "param_best" -name "Umin_best" "toDouble(param_best)"' \
        icmd11='replacecol "param_best" -name "fPDR_best" "toDouble(param_best)"' \
        \
        icmd1='replacecol "param_sigma" -name "Mstar_sigma" "toDouble(param_sigma)"' \
        icmd2='replacecol "param_sigma" -name "LAGN_sigma" "toDouble(param_sigma)"' \
        icmd3='replacecol "param_sigma" -name "LIR_sigma" "toDouble(param_sigma)"' \
        icmd4='replacecol "param_sigma" -name "LIR_warm_sigma" "toDouble(param_sigma)"' \
        icmd5='replacecol "param_sigma" -name "LIR_cold_sigma" "toDouble(param_sigma)"' \
        icmd6='replacecol "param_sigma" -name "Mdust_sigma" "toDouble(param_sigma)"' \
        icmd7='replacecol "param_sigma" -name "Mdust_warm_sigma" "toDouble(param_sigma)"' \
        icmd8='replacecol "param_sigma" -name "Mdust_cold_sigma" "toDouble(param_sigma)"' \
        icmd9='replacecol "param_sigma" -name "Umean_sigma" "toDouble(param_sigma)"' \
        icmd10='replacecol "param_sigma" -name "Umin_sigma" "toDouble(param_sigma)"' \
        icmd11='replacecol "param_sigma" -name "fPDR_sigma" "toDouble(param_sigma)"' \
        \
        icmd1='replacecol "param_L68" -name "Mstar_L68" "toDouble(param_L68)"' \
        icmd2='replacecol "param_L68" -name "LAGN_L68" "toDouble(param_L68)"' \
        icmd3='replacecol "param_L68" -name "LIR_L68" "toDouble(param_L68)"' \
        icmd4='replacecol "param_L68" -name "LIR_warm_L68" "toDouble(param_L68)"' \
        icmd5='replacecol "param_L68" -name "LIR_cold_L68" "toDouble(param_L68)"' \
        icmd6='replacecol "param_L68" -name "Mdust_L68" "toDouble(param_L68)"' \
        icmd7='replacecol "param_L68" -name "Mdust_warm_L68" "toDouble(param_L68)"' \
        icmd8='replacecol "param_L68" -name "Mdust_cold_L68" "toDouble(param_L68)"' \
        icmd9='replacecol "param_L68" -name "Umean_L68" "toDouble(param_L68)"' \
        icmd10='replacecol "param_L68" -name "Umin_L68" "toDouble(param_L68)"' \
        icmd11='replacecol "param_L68" -name "fPDR_L68" "toDouble(param_L68)"' \
        \
        icmd1='replacecol "param_H68" -name "Mstar_H68" "toDouble(param_H68)"' \
        icmd2='replacecol "param_H68" -name "LAGN_H68" "toDouble(param_H68)"' \
        icmd3='replacecol "param_H68" -name "LIR_H68" "toDouble(param_H68)"' \
        icmd4='replacecol "param_H68" -name "LIR_warm_H68" "toDouble(param_H68)"' \
        icmd5='replacecol "param_H68" -name "LIR_cold_H68" "toDouble(param_H68)"' \
        icmd6='replacecol "param_H68" -name "Mdust_H68" "toDouble(param_H68)"' \
        icmd7='replacecol "param_H68" -name "Mdust_warm_H68" "toDouble(param_H68)"' \
        icmd8='replacecol "param_H68" -name "Mdust_cold_H68" "toDouble(param_H68)"' \
        icmd9='replacecol "param_H68" -name "Umean_H68" "toDouble(param_H68)"' \
        icmd10='replacecol "param_H68" -name "Umin_H68" "toDouble(param_H68)"' \
        icmd11='replacecol "param_H68" -name "fPDR_H68" "toDouble(param_H68)"' \
        \
        ocmd='delcols "Source_*"' \
        out="output_datatable_all_in_one.fits" ofmt=fits

if [[ ! -f "output_datatable_all_in_one.fits" ]]; then
    echo "Error occured!"
    exit 255
fi



# 
# Done
# 
echo "Output to \"${output_dir}/output_datatable_all_in_one\"!"
#open "${output_file}"



