#!/bin/bash
# 

set -e


# 
# Define output directory
# 
dir_of_output="SED_fitting_michi2"

if [[ ! -d "$dir_of_output" ]]; then
    echo "Error! Directory \"$dir_of_output\" was not found!"; exit 1
fi

echo cd "$dir_of_output"
cd "$dir_of_output"




# 
# Run python code to limit the SED fitting S/N no higher than the photometry S/N
# 
../$(basename "${BASH_SOURCE[0]}" | sed -e 's/\.sh$//g').py









