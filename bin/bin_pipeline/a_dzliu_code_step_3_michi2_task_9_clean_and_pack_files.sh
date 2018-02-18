#!/bin/bash
# 


cd "SED_fitting_michi2"


find . -name "*.bak" -exec rm \{\} \;

find . -name "*.backup" -exec rm \{\} \;

7z a -t7z -mx9 -ms SED_fitted_files.7z ID_* Results list_*

7z a -t7z -mx9 -ms Results.7z Results list_*







