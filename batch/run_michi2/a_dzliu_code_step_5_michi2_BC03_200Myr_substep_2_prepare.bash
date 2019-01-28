#!/bin/bash
# 

# 
# Define SED_Fitting_Type
# 
SED_Fitting_Type="SED_fitting_michi2_BC03_200Myr"


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

set -e


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
# Prepare parallel running scripts (slurm)
# 
each_task_mem="2G"
each_task_cpu="4"
echo "#!/bin/bash" > "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "#SBATCH --mail-user=dzliu@mpia-hd.mpg.de" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "#SBATCH --ntasks=4 # total task number for sbatch job creation" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "#SBATCH --ntasks-per-core=1 # simultaneous task number per CPU, for sbatch job creation but not srun step creation" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "#SBATCH --cpus-per-task=${each_task_cpu} # multithread" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
if [[ $(hostname) == *"isaac"* ]]; then
echo "#SBATCH --time=4:00:00 # max allowed on isaac is 2-00:00:00" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
elif [[ $(hostname) == *"draco"* ]]; then
echo "#SBATCH --time=4:00:00 # max allowed on draco \"small\" or \"general\" partition is 1-00:00:00" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "#SBATCH --partition=short,small,general # max allowed on draco \"small\" or \"general\" partition is 1-00:00:00" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
else
echo "#SBATCH --time=4:00:00 # try time limit 2-00:00:00" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
fi
echo "#SBATCH --job-name=sedmichi2" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "#SBATCH --output=log_TASK_ID_%a_JOB_ID_%A.out" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
echo "" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"


echo "#!/bin/bash" > "batch_run_${SED_Fitting_Type}_for_all.sh"
echo "#" >> "batch_run_${SED_Fitting_Type}_for_all.sh"


# 
# Loop sources
# 
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
    
    # make working subsubfolder
    if [[ ! -d "${source_name}/${SED_Fitting_Type}" ]]; then
        mkdir -p "${source_name}/${SED_Fitting_Type}"
    fi
    
    # copy photometry data
    if [[ ! -f "${source_name}/${SED_Fitting_Type}/datatable_photometry.txt" ]]; then
        cp "${source_name}/datatable_photometry.txt" \
           "${source_name}/${SED_Fitting_Type}/datatable_photometry.txt"
    fi
    if [[ ! -f "${source_name}/${SED_Fitting_Type}/datatable_id_ra_dec_zprior.txt" ]]; then
        cp "${source_name}/datatable_id_ra_dec_zprior.txt" \
           "${source_name}/${SED_Fitting_Type}/datatable_id_ra_dec_zprior.txt"
    fi
    
    
    echo ""
    echo "Processing ${source_name} ($((i+1))/${#list_of_source_names[@]}) (${SED_Fitting_Type})"
    
    echo "#!/bin/bash" > "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "#" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "if [[ \$(hostname) == \"isaac\"* ]] || [[ \$(hostname) == \"draco\"* ]]; then" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    # load python" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    module load anaconda" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    # also create virtual display" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    if [[ \$(ps aux | grep 'Xvfb :1' | grep -v 'grep' | wc -l) -eq 0 ]]; then" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "        Xvfb :1 -screen 0 1152x900x8 &" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    fi" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    MY_VDISP_PID=\$!" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    export DISPLAY=\":1.0\"" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "fi" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "source \"\$HOME/Cloud/Github/Crab.Toolkit.michi2/SETUP.bash\"" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "current_dir=\$(pwd)" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "echo cd \"${source_name}/${SED_Fitting_Type}\"" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "cd \"${source_name}/${SED_Fitting_Type}\"" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "source_redshift=\$(cat datatable_id_ra_dec_zprior.txt | grep -v '^#' | head -n 1 | sed -e 's/^ *//g' | tr -s ' ' | cut -d ' ' -f 4)" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "michi2-run-SED-fitting-v5 \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    -redshift \$source_redshift \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    -flux datatable_photometry.txt \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    -parallel ${each_task_cpu} \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    -obj-name ${source_name} \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    -lib-stellar BC03.200Myr \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "    -sampling 6000 \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    if [[ ${source_name} == "ID_85001929" ]]; then
    echo "    -no-radio \\" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    else
    echo "    \$@" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    fi
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "echo cd \"\${current_dir}\"" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    echo "cd \"\${current_dir}\"" >> "run_${SED_Fitting_Type}_for_${source_name}.sh"
    # 
    chmod +x "run_${SED_Fitting_Type}_for_${source_name}.sh"
    # 
    echo "srun --exclusive --hint=multithread --unbuffered --cpus-per-task=${each_task_cpu} --ntasks=1 --mem=${each_task_mem} --output=./run_${SED_Fitting_Type}_for_${source_name}.log --error=./run_${SED_Fitting_Type}_for_${source_name}.err ./run_${SED_Fitting_Type}_for_${source_name}.sh &" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"
    # 
    echo "./run_${SED_Fitting_Type}_for_${source_name}.sh" >> "batch_run_${SED_Fitting_Type}_for_all.sh"
    echo "sleep 2.0" >> "batch_run_${SED_Fitting_Type}_for_all.sh"
    # 
done


# 
# write into parallel running scripts (slurm)
# 
echo "wait" >> "batch_slurm_${SED_Fitting_Type}_for_all.sh"


chmod +x "batch_slurm_${SED_Fitting_Type}_for_all.sh"

chmod +x "batch_run_${SED_Fitting_Type}_for_all.sh"
















