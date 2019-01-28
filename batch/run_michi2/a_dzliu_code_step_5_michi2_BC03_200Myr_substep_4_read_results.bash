#!/bin/bash
# 

# 
# Define SED_Fitting_Type
# 
SED_Fitting_Type="SED_fitting_michi2_BC03_200Myr"
param_names=("LIR_total" "Mdust_total" "Mstar" "Umean_total" "LIR_total_40_400")
fit_name="fit_5"
#param_name="LIR_cold"
#param_name="LIR_total_40_400"


# 
# Check software dependencies
# 
if [[ ! -f "$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash" ]]; then
    echo "Error! \"$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash\" was not found! Please clone from \"https://github.com/1054/Crab.Toolkit.michi2\"!"
    exit
fi
source "$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash"


# 
# cd working dir
# 
set -e
echo cd "Multi-wavelength_SEDs"
cd "Multi-wavelength_SEDs"
set +e


# 
# Check necessary files
# 
file_of_source_names="list_of_source_names.txt"
#file_of_source_redshifts="list_of_source_redshifts.txt"
if [[ ! -f "$file_of_source_names" ]]; then echo "Error! \"$file_of_source_names\" was not found!"; exit; fi
#if [[ ! -f "$file_of_source_redshifts" ]]; then echo "Error! \"$file_of_source_redshifts\" was not found!"; exit; fi
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | tr -s ' ' | cut -d ' ' -f 1))
#list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | tr -s ' ' | cut -d ' ' -f 1))


# 
# Make output dir
# 
dir_of_output="../Multi-wavelength_SED_Results/$SED_Fitting_Type"
if [[ "$dir_of_output" != "./" ]]; then
    if [[ ! -d "$dir_of_output" ]]; then
        mkdir -p "$dir_of_output"
    fi
fi
dir_of_output=$(perl -MCwd -e 'print Cwd::abs_path shift' "$dir_of_output")


# 
# Loop sources
# 
icount=0
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    
    source_name="${list_of_source_names[i]}"
    
    # check working subfolder
    if [[ ! -d "${source_name}" ]]; then
        echo "Error! \"$(pwd)/${source_name}\" was not found! Please run substep 1 script first!"
        exit -1
    fi
    
    # check photometry data
    if [[ ! -f "${source_name}/datatable_photometry.txt" ]]; then
        echo "Error! \"$(pwd)/${source_name}/datatable_photometry.txt\" was not found! Please run substep 2 script first!"
        exit -2
    fi
    
    # read redshift
    source_redshift=$(cat "${source_name}/datatable_id_ra_dec_zprior.txt" | grep -v '^#' | head -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 4)
    
    # check working subsubfolder
    if [[ ! -d "${source_name}/${SED_Fitting_Type}" ]]; then
        echo "Error! \"$(pwd)/${source_name}/${SED_Fitting_Type}\" was not found! Please run substep 3 script first!"
        exit -1
    fi
    results_dir="${source_name}/${SED_Fitting_Type}/results_${fit_name}"
	
    
    for param_name in ${param_names[@]}; do 
        # check results for param_name
        if [[ ! -f "${results_dir}/best-fit_param_${param_name}.txt" ]]; then
            echo "Error! \"${results_dir}/best-fit_param_${param_name}.txt\" was not found! Please run substep 3 script first!"
            exit 1
        fi
        # read results for param_name and store into ${dir_of_output}
        if [[ $icount -eq 0 ]]; then
            str_Source=$(printf "# %-10s" "Source")
            str_z=$(printf "%-10s" "z")
            head -n 1 "${results_dir}/best-fit_param_${param_name}.txt" | perl -p -e "s/^/${str_Source} ${str_z} /g" > "${dir_of_output}/best-fit_param_${param_name}.txt"
        fi
        str_Source=$(printf "%-12s" "${source_name}")
        str_z=$(printf "%-10s" "${source_redshift}")
        cat "${results_dir}/best-fit_param_${param_name}.txt" | grep -v '^#' | tail -n 1 | perl -p -e "s/^/${str_Source} ${str_z} /g" >> "${dir_of_output}/best-fit_param_${param_name}.txt"
    done
    
    # copy figures
    cp "${results_dir}/${fit_name}.pdf" "${dir_of_output}/Plot_SED_michi2_${source_name}.pdf"
    cp "${results_dir}/${fit_name}.best.pdf" "${dir_of_output}/Plot_SED_michi2_${source_name}.best.pdf"
    cp "${results_dir}/${fit_name}.chisq.pdf" "${dir_of_output}/Plot_SED_michi2_${source_name}.chisq.pdf"
    cp "${results_dir}/best-fit_SED_${source_name}.txt" "${dir_of_output}/best-fit_SED_${source_name}.txt"
    
    # 
    icount=$((icount+1))
done


echo "Output to \"${dir_of_output}/best-fit_SED_*.txt\"!"
echo "Output to \"${dir_of_output}/best-fit_param_*.txt\"!"
echo "Output to \"${dir_of_output}/Plot_SED_michi2_*.pdf\"!"






