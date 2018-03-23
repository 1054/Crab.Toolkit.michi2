#!/bin/bash
# 

set -e


# 
# Define output directory
# 
dir_of_output="."

if [[ ! -d "$dir_of_output" ]]; then
    echo "Error! Directory \"$dir_of_output\" was not found!"; exit 1
fi

echo cd "$dir_of_output"
cd "$dir_of_output"




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
list_of_source_names=($(cat "$file_of_source_names" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'| tr -s ' ' | cut -d ' ' -f 1))
list_of_source_redshifts=($(cat "$file_of_source_redshifts" | grep -v '^#' | sed -e 's/^ *//g' | sed -e 's/ *$//g'| tr -s ' ' | cut -d ' ' -f 1))


# 
# Define second-level output directory
# 
output_dir="Results"


# 
# make output directory and write ReadMe file
# 
if [[ ! -d $output_dir ]]; then mkdir $output_dir; fi

echo "The files in this directory (or zip package) are the results of multi-component SED fitting." >  $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "ID_*_fit_5_components.pdf           -- This figure shows a range of fitted SEDs, " >> $output_dir/ReadMe
echo "                                       each single SED component and the observed photometry data points." >> $output_dir/ReadMe
echo "ID_*_fit_5_components.chisq.pdf     -- This figure shows the chi-square distributions of fitted parameters. " >> $output_dir/ReadMe
echo "                                       Each parameter has two panels, first with larger plotting ranges and second zoomed in." >> $output_dir/ReadMe
echo "ID_*_fit_5_component_SEDs/SED_LIB*  -- These are the best-fit single component SEDs." >> $output_dir/ReadMe
echo "                                       (Note that X is rest-frame and needs to be multiplied by (1+z) for displaying, " >> $output_dir/ReadMe
echo "                                       while no (1+z) is needed for Y. See ReadMe therein.)" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "The 5 components we fitted are: " >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "    LIB1     -- BC03 (200Myr, solar metallicity, constant SFH, Chabrier IMF) SEDs with Calzetti+2000 attenuation law. Free parameter: E(B-V)." >> $output_dir/ReadMe
echo "    LIB2     -- Mid-infrared AGN SEDs based on observations, from Mullaney+2011. Free parameter: AGN_Type." >> $output_dir/ReadMe
echo "    LIB3     -- DL07 warm dust SEDs. Free parameters: Umin, qPAH." >> $output_dir/ReadMe
echo "    LIB4     -- DL07 cold dust SEDs. Free parameters: Umin, qPAH. (Note that the Umin and qPAH are locked to be the same as in LIB3 while fitting.)" >> $output_dir/ReadMe
echo "    LIB5     -- Radio SED, simple power-law with slope -0.8, fixed to L_IR(8-1000um,LIB3+LIB4) via the IR-radio correlation with qIR=2.4. Currently no parameter." >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "The parameters in the chi-square distribution figure are: " >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "    log M_*                         -- Stellar mass, in unit of solar mass." >> $output_dir/ReadMe
echo "    log L_AGN                       -- AGN luminosity (TODO: needs further confirmation on the fitted absolute luminosities)" >> $output_dir/ReadMe
echo "    U_min (cold, warm)              -- The U_min parameter in DL07 dust model, which is the ambient interstellar radiation field strength." >> $output_dir/ReadMe
echo "    log L_IR (cold, warm, total)    -- The L_IR(8-1000um) luminosity from cold/warm/cold+warm dust." >> $output_dir/ReadMe
echo "    log M_dust (cold, warm, total)  -- The mass of cold/warm/cold+warm dust, in unit of solar mass." >> $output_dir/ReadMe
echo "    log \delta_{PDR} (total)        -- The fraction of warm dust (coming from photon-dominated regions, PDRs) mass to the total (cold+warm) dust mass." >> $output_dir/ReadMe
echo "    <U> (total)                     -- The mean interstellar radiation field (ambient+PDR) strength of the whole galaxy." >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "We also output the best-fit data table for several parameters: " >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "    best-fit_param_*.txt" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "Note that the columns in each table are: " >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "    Source        -- The source name" >> $output_dir/ReadMe
echo "    *_median      -- The median of high probability parameter values" >> $output_dir/ReadMe
echo "    *_best        -- The minimum chi-square solution's parameter value" >> $output_dir/ReadMe
echo "    *_sigma       -- The 1-sigma uncertainty of the fitted parameter value" >> $output_dir/ReadMe
echo "    *_L68         -- The lower 68% confidence parameter value" >> $output_dir/ReadMe
echo "    *_H68         -- The upper 68% confidence parameter value" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
echo "" >> $output_dir/ReadMe
# \"michi2\" SED fitting (https://github.com/1054/Crab.Toolkit.michi2; Liu et al. 2018 in prep.)

# 
# first determine the max string length of the source name
# 
_SourceNX=8 # should at least longer than "# Source" string.
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    _Source="${list_of_source_names[i]}"
    if [[ ${#_Source} -gt $_SourceNX ]]; then _SourceNX=${#_Source}; fi
done


# 
# then copy best-fit figures and combine data tables
# 
for (( i = 0; i < ${#list_of_source_names[@]}; i++ )); do
    
    _Source="${list_of_source_names[i]}"
    
    echo "$_Source"
    
    # 
    # Copy SED fitting figures
    cp ./$_Source/fit_5.pdf           $output_dir/${_Source}_fit_5_components.pdf
    cp ./$_Source/fit_5.best.pdf      $output_dir/${_Source}_fit_5_components.best.pdf
    cp ./$_Source/fit_5.chisq.pdf     $output_dir/${_Source}_fit_5_components.chisq.pdf
    
    if [[ -f ./$_Source/fit_2.pdf ]]; then cp ./$_Source/fit_2.pdf $output_dir/${_Source}_fit_2_components.pdf; fi
    if [[ -f ./$_Source/fit_2.chisq.pdf ]]; then cp ./$_Source/fit_2.pdf $output_dir/${_Source}_fit_2_components.chisq.pdf; fi
    
    # 
    # Copy SED fitting best-fit ascii file
    if [[ ! -d $output_dir/${_Source}_fit_5_component_SEDs ]]; then 
        mkdir -p $output_dir/${_Source}_fit_5_component_SEDs
    fi
    cp ./$_Source/obj_1/SED_LIB[12345]   $output_dir/${_Source}_fit_5_component_SEDs/
    cp ./$_Source/obj_1/SED_SUM          $output_dir/${_Source}_fit_5_component_SEDs/
    echo "X is the rest-frame wavelength in unit of um."                     >  $output_dir/${_Source}_fit_5_component_SEDs/ReadMe
    echo "Y is the flux density in unit of mJy."                             >> $output_dir/${_Source}_fit_5_component_SEDs/ReadMe
    echo "Redshift (1+z) item has not been multiplied to the X array!"       >> $output_dir/${_Source}_fit_5_component_SEDs/ReadMe
    echo "To plot the SED in obs-frame, just do like \"plot(X*(1+z),Y)\"."   >> $output_dir/${_Source}_fit_5_component_SEDs/ReadMe
    
    # 
    # Combine best-fit parameter data tables
    for _Param in LIR_total Mdust_total Umean_total LAGN Mstar fPDR_total; do
        # Dealing with column header string length
        _ParN1=$(awk "BEGIN { if (${#_Param}>=5) print ${#_Param}-5; else print 0; }")
        _ParN2=$(awk "BEGIN { if (${#_Param}>=5) print 0; else print 5-${#_Param}; }")
        if [[ $_ParN1 -gt 0 ]]; then _ParW1="$(seq -s ' ' $_ParN1 | tr -d '[:digit:]')"; else _ParW1=""; fi # padding white space
        if [[ $_ParN2 -gt 0 ]]; then _ParW2="$(seq -s ' ' $_ParN2 | tr -d '[:digit:]')"; else _ParW2=""; fi # padding white space
        #echo "${_ParW1}param_"
        #echo "${_ParW2}${_Param}_"
        # 
        # Now write to file
        _SourceN1=$(awk "BEGIN {if (8>${_SourceNX}) print 0; else print (${_SourceNX}-8);}") # count how many white space to pad for the string "# Source"
        _SourceN2=$(awk "BEGIN {if (${#_Source}>${_SourceNX}) print 0; else print (${_SourceNX}-${#_Source});}") # count how many white space to pad for the string "${_Source}"
        if [[ $_SourceN1 -gt 0 ]]; then _SourceW1=$(seq -s ' ' $_SourceN1 | tr -d '[:digit:]'); else _SourceW1=""; fi # padding white space
        if [[ $_SourceN2 -gt 0 ]]; then _SourceW2=$(seq -s ' ' $_SourceN2 | tr -d '[:digit:]'); else _SourceW2=""; fi # padding white space
        if [[ $i -eq 0 ]]; then
        head -n 1 ./$_Source/best-fit_param_$_Param.txt  | sed -e "s/^/# Source${_SourceW1} /g" | sed -e "s/${_ParW1}param_/${_ParW2}${_Param}_/g"  >   $output_dir/best-fit_param_$_Param.txt
        fi
        tail -n 1 ./$_Source/best-fit_param_$_Param.txt  | sed -e "s/^/$_Source${_SourceW2} /g"                                                     >>  $output_dir/best-fit_param_$_Param.txt
        #echo "--# Source${_SourceW1}--"
        #echo "--$_Source${_SourceW2}--"
    done
    
    
    # 
    # get chisq <20180216>
    if [[ $i -eq 0 ]]; then
    echo "# Source${_SourceW1}  chisq" > $output_dir/best-fit_chisq.txt
    fi
    echo "$_Source${_SourceW2}  "$(cat ./$_Source/obj_1/chi2.txt | grep -v '^#' | head -n 1) >> $output_dir/best-fit_chisq.txt
    
    
done





