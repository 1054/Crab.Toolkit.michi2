#!/bin/bash
# 

echo "macro read makelibSED_DL07_spec_2p0.sm makelibSED_DL07" | sm | tee makelibSED_DL07_spec_2p0.log

7z a -tzip lib.DL07UPD2010.DuoCom.SED.zip lib.DL07.2010.03.18.spec.2.0.HiExCom.SED lib.DL07.2010.03.18.spec.2.0.LoExCom.SED

#rm DL07spec

