#!/bin/bash
#SBATCH --mail-user=dzliu@mpia-hd.mpg.de
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --time=24:00:00
#SBATCH --mem=4000
#SBATCH --cpus-per-task=32
#SBATCH --output=log_SED_fitting_michi2_TASK_ID_%a_JOB_ID_%A.out


# 
# What does this script do --
# This script runs michi2 SED fitting
# 


# 
# To run this script on the ISAAC machine -- 
# srun -N1 -n10 -l a_dzliu_code_for_running_SED_fitting_michi2_on_ISAAC.sh
# or
# sbatch a_dzliu_code_for_running_SED_fitting_michi2_on_ISAAC.sh
# or 
# srun -N1 -n20 --pty bash
# or
# sbatch --array=1-6%2 a_dzliu_code_step_2_run_michi2_task_3_run_sbatch_on_isaac.sh
# 


# 
# to run this script in Slurm job array mode
# sbatch --array=1-1%1 -N1 ~/Cloud/Github/AlmaCosmos/Pipeline/a3cosmos-MC-simulation-physically-motivated/do_simulation/a_dzliu_code_for_Simulation_on_ISAAC_Step_2_Simulate.sh test
# sbatch --array=4-4%1 -N1 ~/Cloud/Github/AlmaCosmos/Pipeline/a3cosmos-MC-simulation-physically-motivated/do_simulation/a_dzliu_code_for_Simulation_on_ISAAC_Step_2_Simulate.sh
# 
echo "Hostname: "$(/bin/hostname)
echo "PWD: "$(/bin/pwd)
#echo "SLURM_JOBID: "$SLURM_JOBID
#echo "SLURM_JOB_NODELIST: "$SLURM_JOB_NODELIST
#echo "SLURM_NNODES: "$SLURM_NNODES
#echo "SLURM_ARRAY_TASK_ID: "$SLURM_ARRAY_TASK_ID
#echo "SLURM_ARRAY_JOB_ID: "$SLURM_ARRAY_JOB_ID
#echo "SLURMTMPDIR: "$SLURMTMPDIR
#echo "SLURM_SUBMIT_DIR: "$SLURM_SUBMIT_DIR


# 
# Check user input
# 
if [[ $# -eq 0 ]]; then
    echo "Usage: "
    echo "    Please input source name! E.g., \"ID_1234\""
    echo ""
    exit
fi


# 
# Load module
# 
if [[ $(hostname) == isaac* ]]; then
    module load anaconda
fi


# 
# Check software dependencies
# 
Github_dir="$HOME/Cloud/Github"
if [[ $(hostname) == isaac* ]]; then
    Github_dir="/u/$USER/Cloud/Github"
fi

#if [[ ! -f "$Github_dir/DeepFields.SuperDeblending/Softwares/SETUP.bash" ]]; then
#    echo "Error! \"$Github_dir/DeepFields.SuperDeblending/Softwares/SETUP.bash\" was not found! Please clone from \"https://github.com/1054/DeepFields.SuperDeblending\"!"
#    exit
#else
#    source "$Github_dir/DeepFields.SuperDeblending/Softwares/SETUP.bash"
#fi

if [[ ! -f "$Github_dir/Crab.Toolkit.michi2/SETUP.bash" ]]; then
    echo "Error! \"$Github_dir/Crab.Toolkit.michi2/SETUP.bash\" was not found! Please clone from \"https://github.com/1054/Crab.Toolkit.michi2\"!"
    exit
else
    source "$Github_dir/Crab.Toolkit.michi2/SETUP.bash"
fi


# 
# CD
# 

if [[ ! -f "list_of_source_names.txt" ]] || [[ ! -f "list_of_source_redshifts.txt" ]]; then
    echo "Error! \"list_of_source_names.dat\" or \"list_of_source_redshifts.txt\" was not found!"
    exit 1
fi
list_of_source_names=($(cat "list_of_source_names.txt" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))
list_of_source_redshifts=($(cat "list_of_source_redshifts.txt" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'))


# 
# Run SED fitting
# 
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    
    itask=$((i+1))
    
    if [[ x"${list_of_source_names[i]}" == x"" ]]; then
        continue
    fi
    
    # check user input, only run for sources input by the user.
    imatch=0
    for iargs in "$@"; do
        if [[ "${list_of_source_names[i]}" == $(echo "${iargs}" | sed -e 's%/$%%g') ]]; then
            imatch=1
            break
        fi
    done
    if [[ imatch -eq 0 ]]; then
        continue
    fi
    
    
    if [[ ! -d "${list_of_source_names[i]}" ]]; then
        echo "Warning! \"${list_of_source_names[i]}\" was not found! Will skip this source!"
        continue
    fi
    
    cd "${list_of_source_names[i]}"
    
    
    if [[ ! -f "fit_5.out" ]]; then
        
        #michi2-deploy-files
        
        michi2-run-fitting-5-components ${list_of_source_redshifts[i]} -flux extracted_flux.txt -parallel 16
        
        rm -rf obj_* best* Plot_*
        
        michi2-plot-fitting-results fit_5.out -flux extracted_flux.txt -source "${list_of_source_names[i]}"
        
        sleep 3
    
    else
        
        echo "Found existing \"fit_5.out\" under directory \"${list_of_source_names[i]}\"! Will skip this source!"
    
    fi
    
    cd "../"
    
done









