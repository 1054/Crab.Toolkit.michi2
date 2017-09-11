#!/bin/bash
#


if [[ ! -f "bestfit_Umean_range.txt" ]]; then
    echo "Error! \"bestfit_Umean_range.txt\" was not found!"; exit
fi
if [[ ! -f "bestfit_EBV_range.txt" ]]; then
    echo "Error! \"bestfit_EBV_range.txt\" was not found!"; exit
fi
if [[ ! -f "bestfit_LTIR_range_log10.txt" ]]; then
    echo "Error! \"bestfit_LTIR_range_log10.txt\" was not found!"; exit
fi
if [[ ! -f "bestfit_Mdust_range_log10.txt" ]]; then
    echo "Error! \"bestfit_Mdust_range_log10.txt\" was not found!"; exit
fi
if [[ ! -f "bestfit_Mstar_range_log10.txt" ]]; then
    echo "Error! \"bestfit_Mstar_range_log10.txt\" was not found!"; exit
fi


echo "macro read run_printing_bestfit.sm print_bestfit" | sm

