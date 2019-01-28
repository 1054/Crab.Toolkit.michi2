#!/bin/bash
# 

# 
# Define SED_Fitting_Type
# 
SED_Fitting_Type="SED_fitting_michi2_BC03_200Myr"



# 
# Submit the slurm job
#
cd "Multi-wavelength_SEDs/"
sbatch "batch_slurm_${SED_Fitting_Type}_for_all.sh"

















