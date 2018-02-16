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





mkdir "Results"


for (( i=0; i<${#list_of_source_names[@]}; i++ )); do
    
    if [[ 1 == 1 ]]; then
        # 
        echo 
        echo 
        echo "***********************************************"
        echo "cd \"${list_of_source_names[i]}/\""
        echo "***********************************************"
        cd "${list_of_source_names[i]}/"
        source_id=$(echo "${list_of_source_names[i]}" | sed -e 's/ID_//g')
        # 
        if [[ $i -eq 0 ]]; then
        printf "# %-20s %15s\n" "Source" "z" | sed -e "s/$/$(head -n 1 fit_parallel_HDFN/ResLMT_${source_id}.txt | sed -e 's/^#/ /g')/g" > "../Results/best-fit_params.txt"
        fi
        printf "# %-20s %15s\n" "Source" "z" | sed -e "s/$/$(head -n 1 fit_parallel_HDFN/ResLMT_${source_id}.txt | sed -e 's/^#/ /g')/g"
        printf "  %-20s %15.5f\n" "${list_of_source_names[i]}" "${list_of_source_redshifts[i]}" | sed -e "s/$/$(cat fit_parallel_HDFN/ResLMT_${source_id}.txt | grep -v '^#' | head -n 1)/g"
        printf "  %-20s %15.5f\n" "${list_of_source_names[i]}" "${list_of_source_redshifts[i]}" | sed -e "s/$/$(cat fit_parallel_HDFN/ResLMT_${source_id}.txt | grep -v '^#' | head -n 1)/g" >> "../Results/best-fit_params.txt"
        # 
        if [[ $i -eq 0 ]]; then
        printf "# %-20s %15s\n" "Source" "z" | sed -e "s/$/$(head -n 1 fit_parallel_HDFN/fit_${source_id}.csv | sed -e 's/^#/ /g')/g" > "../Results/best-fit_chisq.txt"
        fi
        printf "# %-20s %15s\n" "Source" "z" | sed -e "s/$/$(head -n 1 fit_parallel_HDFN/fit_${source_id}.csv | sed -e 's/^#/ /g')/g"
        printf "  %-20s %15.5f\n" "${list_of_source_names[i]}" "${list_of_source_redshifts[i]}" | sed -e "s/$/$(cat fit_parallel_HDFN/fit_${source_id}.csv | grep -v '^#' | head -n 1)/g"
        printf "  %-20s %15.5f\n" "${list_of_source_names[i]}" "${list_of_source_redshifts[i]}" | sed -e "s/$/$(cat fit_parallel_HDFN/fit_${source_id}.csv | grep -v '^#' | head -n 1)/g" >> "../Results/best-fit_chisq.txt"
        # 
        echo "cd \"../\""
        cd "../"
    fi
    
    #break
    
done


