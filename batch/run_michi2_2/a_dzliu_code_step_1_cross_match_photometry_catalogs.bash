#!/bin/bash
# 

#if [[ ! -f datatable_ID_RA_Dec_CO_detected.txt ]]; then
#    
#fi

if [[ ! -d datatable_photometry ]]; then
    mkdir datatable_photometry
fi
echo cd datatable_photometry
cd datatable_photometry

set -e

input_cat_ID_RA_Dec="../datatable_ID_RA_Dec.txt"
ifmt_cat_ID_RA_Dec=ascii
input_cat_dir="/Catalogs"
input_cat_Laigle2016="$input_cat_dir/COSMOS_photometry_Laigle2016/COSMOS2015_Laigle+_v1.1.fits"
input_cat_Jin2018="$input_cat_dir/COSMOS_photometry_Jin2017_SuperDeblend_FIR_mm/v20180719/COSMOS_Super_Deblended_FIRmm_Catalog_20180719.fits"
input_cat_Ashby2018="$input_cat_dir/A3COSMOS_Master_Catalog_20181106/Catalog_Ashby_2018_SMUVS_all_combined_by_dzliu.fits.gz"
input_cat_Liu2020a="$input_cat_dir/A3COSMOS_20180801/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_SED_fluxes.fits"

ls "$input_cat_Laigle2016"
ls "$input_cat_Jin2018"
ls "$input_cat_Ashby2018"
ls "$input_cat_Liu2020a"

crossmatched_seplimit=1.4 # arcsec

# 
# Convert input cat ID RA Dec
# 
i=0
echo "************"
echo "*** ${i} ***"
echo "************"
if [[ ! -f "datatable_crossmatched_${i}.fits" ]]; then
topcat -stilts tpipe \
                in="$input_cat_ID_RA_Dec" \
                ifmt="$ifmt_cat_ID_RA_Dec" \
                cmd='keepcols "ID RA Dec z"' \
                out="datatable_crossmatched_${i}.fits" \
                ofmt=fits
topcat -stilts tpipe \
                in="datatable_crossmatched_${i}.fits" \
                ifmt=fits \
                omode=meta \
               | tee "datatable_crossmatched_${i}.meta.txt"
               # meta
fi


# 
# Now we cross-match the Laigle2016 Catalog only for photometry data
# 
i=$((i+1))
echo "************"
echo "*** ${i} ***"
echo "************"
if [[ ! -f "datatable_crossmatched_${i}.fits" ]]; then
topcat -stilts tmatchn \
                nin=2 \
                in1="datatable_crossmatched_$((i-1)).fits" \
                ifmt1=fits \
                values1="RA Dec" \
                suffix1="" \
                join1=always \
                in2="$input_cat_Laigle2016" \
                ifmt2=fits \
                icmd2='addcol -before $1 ID NUMBER' \
                icmd2='addcol -after ID RA -desc "corrected for astrometry" "ALPHA_J2000 - 0.09/3600.0/cos(DELTA_J2000/180.0*PI)"' \
                icmd2='addcol -after RA Dec -desc "corrected for astrometry" "DELTA_J2000 - 0.015/3600.0"' \
                icmd2='delcols "*_IMAGE FLUX_RADIUS KRON_RADIUS NAME_* *XMM* *CHANDRA* *_20CM *_90CM *GALEX*"' \
                icmd2='delcols "FLAG_XRAYBLEND OFFSET PHOTOZ TYPE ZPDF ZPDF_L68 ZPDF_H68 ZMINCHI2 CHI2_BEST ZP_2 CHI2_2 NBFILT ZQ CHIQ MODQ MODS CHIS MODEL AGE EXTINCTION MNUV MU MB MV MR MI MZ MY MJ MH MK MNUV_MR CLASS MASS_MED MASS_MED_MIN68 MASS_MED_MAX68 MASS_BEST SFR_MED SFR_MED_MIN68 SFR_MED_MAX68 SFR_BEST SSFR_MED SSFR_MED_MIN68 SSFR_MED_MAX68 SSFR_BEST L_NU L_R L_K"' \
                icmd2='delcols "*_APER2 *_APER2_* *_MAG *_MAG_* *_MAGERR *_MAGERR_* *_FLAGS *_FLAGS_* *_IMAFLAGS *_IMAFLAGS_* *_NUSTAR *_NUSTAR_* ID_A24 ID2006 ID2008 ID2013"' \
                icmd2='delcols "FLUX_814W FLUXERR_814W"' \
                icmd2='delcols "FLUX_24 FLUXERR_24 MAG_24 MAGERR_24 FLUX_100 FLUXERR_100 FLUX_160 FLUXERR_160 FLUX_250 FLUXERR_250 FLUXERRTOT_250 FLUX_350 FLUXERR_350 FLUXERRTOT_350 FLUX_500 FLUXERR_500 FLUXERRTOT_500"' \
                values2="RA Dec" \
                suffix2="_Laigle2016" \
                multimode=pairs \
                matcher=sky \
                params="${crossmatched_seplimit}" \
                fixcols=all \
                ocmd='addcol Sep_Laigle2016 -units "arcsec" "sqrt(pow((RA-RA_Laigle2016)*cos(Dec_Laigle2016/180.0*PI),2)+pow(Dec-Dec_Laigle2016,2))*3600.0"' \
                ofmt=fits \
                out="datatable_crossmatched_${i}.fits"
                #<TODO>#
                #icmd2='delcols "FLUX_814W FLUXERR_814W"' \
                #icmd2='replacecol FLUX_814W -name FLUX_814W -units "uJy" "FLUX_814W > -90 ? pow(10, FLUX_814W/(-2.5) ) * 3630.780548 * 1e6 : -99"' \
                #icmd2='replacecol FLUXERR_814W -name FLUXERR_814W -units "uJy" "FLUX_814W > -90 ? FLUXERR_814W * FLUX_814W : -99"' \
                # 
topcat -stilts tpipe \
                in="datatable_crossmatched_${i}.fits" \
                ifmt=fits \
                omode=meta \
               | tee "datatable_crossmatched_${i}.meta.txt"
               # meta
fi


# 
# Now we cross-match the Jin2018 Catalog only for photometry data
# 
i=$((i+1))
echo "************"
echo "*** ${i} ***"
echo "************"
if [[ ! -f "datatable_crossmatched_${i}.fits" ]]; then
topcat -stilts tmatchn \
                nin=2 \
                in1="datatable_crossmatched_$((i-1)).fits" \
                ifmt1=fits \
                values1="RA Dec" \
                suffix1="" \
                join1=always \
                in2="$input_cat_Jin2018" \
                ifmt2=fits \
                icmd2="replacecol RA -name \"RA\" -desc \"corrected for astrometry for Laigle2016 priors\" \"(ID<1e8) ? RA - 0.09/3600.0/cos(DEC/180.0*PI) : RA\"" \
                icmd2="replacecol DEC -name \"Dec\" -desc \"corrected for astrometry for Laigle2016 priors\" \"(ID<1e8) ? DEC - 0.015/3600.0 : DEC\"" \
                icmd2="replacecol ID -name \"ID\" \"ID\"" \
                icmd2="addcol FLUX_K -before Kmag -units \"mJy\" \"(Kmag!=-99) ? pow(10,(-0.4*(Kmag-23.9)-3.0)) : -99\"" \
                icmd2="addcol FLUXERR_K -before Kmag -units \"mJy\" \"(Kmag!=-99) ? FLUX_K/10.0 : 1e10\"" \
                icmd2="replacecol DF24 -name \"FLUXERR_24\" -units \"mJy\" \"(F24!=-99) ? DF24 : 1e10\"" \
                icmd2="replacecol F24 -name \"FLUX_24\" -units \"mJy\" \"(F24!=-99) ? F24 : -99\"" \
                icmd2="replacecol DFCH1 -name \"FLUXERR_IRAC1\" -units \"mJy\" \"(FCH1!=-99) ? DFCH1 : 1e10\"" \
                icmd2="replacecol FCH1 -name \"FLUX_IRAC1\" -units \"mJy\" \"(FCH1!=-99) ? FCH1 : -99\"" \
                icmd2="replacecol DFCH2 -name \"FLUXERR_IRAC2\" -units \"mJy\" \"(FCH2!=-99) ? DFCH2 : 1e10\"" \
                icmd2="replacecol FCH2 -name \"FLUX_IRAC2\" -units \"mJy\" \"(FCH2!=-99) ? FCH2 : -99\"" \
                icmd2="replacecol DFCH3 -name \"FLUXERR_IRAC3\" -units \"mJy\" \"(FCH3!=-99) ? DFCH3 : 1e10\"" \
                icmd2="replacecol FCH3 -name \"FLUX_IRAC3\" -units \"mJy\" \"(FCH3!=-99) ? FCH3 : -99\"" \
                icmd2="replacecol DFCH4 -name \"FLUXERR_IRAC4\" -units \"mJy\" \"(FCH4!=-99) ? DFCH4 : 1e10\"" \
                icmd2="replacecol FCH4 -name \"FLUX_IRAC4\" -units \"mJy\" \"(FCH4!=-99) ? FCH4 : -99\"" \
                icmd2="replacecol DF10CM -name \"FLUXERR_10cm\" -units \"mJy\" \"(F10CM!=-99) ? DF10CM : 1e10\"" \
                icmd2="replacecol F10CM -name \"FLUX_10cm\" -units \"mJy\" \"(F10CM!=-99) ? F10CM : -99\"" \
                icmd2="replacecol DF20CM -name \"FLUXERR_20cm\" -units \"mJy\" \"(F20CM!=-99) ? DF20CM : 1e10\"" \
                icmd2="replacecol F20CM -name \"FLUX_20cm\" -units \"mJy\" \"(F20CM!=-99) ? F20CM : -99\"" \
                icmd2="replacecol F100 -name \"FLUX_100\" -units \"mJy\" \"F100\"" \
                icmd2="replacecol F160 -name \"FLUX_160\" -units \"mJy\" \"F160\"" \
                icmd2="replacecol F250 -name \"FLUX_250\" -units \"mJy\" \"F250\"" \
                icmd2="replacecol F350 -name \"FLUX_350\" -units \"mJy\" \"F350\"" \
                icmd2="replacecol F500 -name \"FLUX_500\" -units \"mJy\" \"F500\"" \
                icmd2="replacecol F850 -name \"FLUX_850\" -units \"mJy\" \"F850\"" \
                icmd2="replacecol F1100 -name \"FLUX_1100\" -units \"mJy\" \"F1100\"" \
                icmd2="replacecol F1200 -name \"FLUX_1200\" -units \"mJy\" \"F1200\"" \
                icmd2="replacecol DF100 -name \"FLUXERR_100\" -units \"mJy\" \"DF100\"" \
                icmd2="replacecol DF160 -name \"FLUXERR_160\" -units \"mJy\" \"DF160\"" \
                icmd2="replacecol DF250 -name \"FLUXERR_250\" -units \"mJy\" \"DF250\"" \
                icmd2="replacecol DF350 -name \"FLUXERR_350\" -units \"mJy\" \"DF350\"" \
                icmd2="replacecol DF500 -name \"FLUXERR_500\" -units \"mJy\" \"DF500\"" \
                icmd2="replacecol DF850 -name \"FLUXERR_850\" -units \"mJy\" \"DF850\"" \
                icmd2="replacecol DF1100 -name \"FLUXERR_1100\" -units \"mJy\" \"DF1100\"" \
                icmd2="replacecol DF1200 -name \"FLUXERR_1200\" -units \"mJy\" \"DF1200\"" \
                icmd2='replacecol FLUXERR_IRAC1 "(FLUXERR_IRAC1>=1e6) ? 1e10 : FLUXERR_IRAC1"' \
                icmd2='replacecol FLUXERR_IRAC2 "(FLUXERR_IRAC2>=1e6) ? 1e10 : FLUXERR_IRAC2"' \
                icmd2='replacecol FLUXERR_IRAC3 "(FLUXERR_IRAC3>=1e6) ? 1e10 : FLUXERR_IRAC3"' \
                icmd2='replacecol FLUXERR_IRAC4 "(FLUXERR_IRAC4>=1e6) ? 1e10 : FLUXERR_IRAC4"' \
                icmd2='replacecol FLUXERR_K "(FLUXERR_K>=1e3) || (FLUX_K<1e-6) ? 1e10 : FLUXERR_K"' \
                icmd2='replacecol FLUXERR_24 "(FLUXERR_24>=1e3) || (FLUX_24<1e-6) ? 1e10 : FLUXERR_24"' \
                icmd2='replacecol FLUXERR_20cm "(FLUXERR_20cm>=1e3) || (FLUX_20cm<1e-6) ? 1e10 : FLUXERR_20cm"' \
                icmd2='replacecol FLUXERR_100 "(FLUX_100<3.0*FLUXERR_100 && FLUXERR_100>=10) || (FLUX_100<1e-4) ? 1e10 : FLUXERR_100"' \
                icmd2='replacecol FLUXERR_160 "(FLUX_160<3.0*FLUXERR_160 && FLUXERR_160>=15) || (FLUX_160<1e-4) ? 1e10 : FLUXERR_160"' \
                icmd2='replacecol FLUXERR_250 "(FLUX_250<3.0*FLUXERR_250 && FLUXERR_250>=15) || (FLUX_250<1e-5) ? 1e10 : FLUXERR_250"' \
                icmd2='replacecol FLUXERR_350 "(FLUX_350<3.0*FLUXERR_350 && FLUXERR_350>=20) || (FLUX_350<1e-5) ? 1e10 : FLUXERR_350"' \
                icmd2='replacecol FLUXERR_500 "(FLUX_500<3.0*FLUXERR_500 && FLUXERR_500>=20) || (FLUX_500<1e-5) ? 1e10 : FLUXERR_500"' \
                icmd2='replacecol FLUXERR_850 "(FLUX_850<3.0*FLUXERR_850 && FLUXERR_850>=10) || (FLUX_850<1e-5) ? 1e10 : FLUXERR_850"' \
                icmd2='replacecol FLUXERR_1100 "(FLUX_1100<3.0*FLUXERR_1100 && FLUXERR_1100>=10) || (FLUX_1100<1e-5) ? 1e10 : FLUXERR_1100"' \
                icmd2='replacecol FLUXERR_1200 "(FLUX_1200<3.0*FLUXERR_1200 && FLUXERR_1200>=10) || (FLUX_1200<1e-5) ? 1e10 : FLUXERR_1200"' \
                icmd2="delcols \"Kmag F10CM_* DF10CM_* \"" \
                values2="RA Dec" \
                suffix2="_Jin2018" \
                multimode=pairs \
                matcher=sky \
                params="${crossmatched_seplimit}" \
                fixcols=all \
                ocmd='addcol Sep_Jin2018 -units "arcsec" "sqrt(pow((RA-RA_Jin2018)*cos(Dec_Jin2018/180.0*PI),2)+pow(Dec-Dec_Jin2018,2))*3600.0"' \
                ofmt=fits \
                out="datatable_crossmatched_${i}.fits"
                # 
topcat -stilts tpipe \
                in="datatable_crossmatched_${i}.fits" \
                ifmt=fits \
                omode=meta \
               | tee "datatable_crossmatched_${i}.meta.txt"
               # meta
fi


# 
# Now we cross-match the Ashby2018 Catalog only for photometry data
# here we need backward cross-match because the IRAC image is fitted with PSF!
# 
i=$((i+1))
echo "************"
echo "*** ${i} ***"
echo "************"
if [[ ! -f "datatable_crossmatched_${i}.fits" ]]; then
topcat -stilts tmatchn \
                nin=2 \
                in1="datatable_crossmatched_$((i-1)).fits" \
                ifmt1=fits \
                values1="RA Dec" \
                suffix1="" \
                join1=always \
                in2="$input_cat_Ashby2018" \
                icmd2='addcol FLUX_IRAC1 -units "mJy" "(MAG1_PSF>0 && MAG1_PSF<99) ? pow(10,MAG1_PSF/(-2.5)) * 3630.780548 * 1e3 : -99"' \
                icmd2='addcol FLUX_IRAC2 -units "mJy" "(MAG2_PSF>0 && MAG2_PSF<99) ? pow(10,MAG2_PSF/(-2.5)) * 3630.780548 * 1e3 : -99"' \
                icmd2='addcol FLUXERR_IRAC1 -units "mJy" "(MAG1_PSF>0 && MAG1_PSF<99) ? MAG1_ERR * pow(10,MAG1_PSF/(-2.5)) * 3630.780548 * 1e3 : 1e10"' \
                icmd2='addcol FLUXERR_IRAC2 -units "mJy" "(MAG2_PSF>0 && MAG2_PSF<99) ? MAG2_ERR * pow(10,MAG2_PSF/(-2.5)) * 3630.780548 * 1e3 : 1e10"' \
                icmd2='addcol ID "INDEX"' \
                icmd2='replacecol ra -name "RA" "ra"' \
                icmd2='replacecol dec -name "Dec" "dec"' \
                icmd2='replacecol Flag1 -name "Flag_IRAC1" -desc "1 means possibly contaminated by bright star" Flag1' \
                icmd2='replacecol Flag2 -name "Flag_IRAC2" -desc "1 means possibly contaminated by bright star" Flag2' \
                icmd2='keepcols "ID RA Dec FLUX_IRAC1 FLUXERR_IRAC1 Flag_IRAC1 FLUX_IRAC2 FLUXERR_IRAC2 Flag_IRAC2"' \
                values2="RA Dec" \
                suffix2="_Ashby2018" \
                multimode=pairs \
                matcher=sky \
                params="${crossmatched_seplimit}" \
                fixcols=all \
                ocmd='addcol Sep_Ashby2018 -units "arcsec" "sqrt(pow((RA-RA_Ashby2018)*cos(Dec_Ashby2018/180.0*PI),2)+pow(Dec-Dec_Ashby2018,2))*3600.0"' \
                ofmt=fits \
                out="datatable_crossmatched_${i}.fits"
                # 
topcat -stilts tpipe \
                in="datatable_crossmatched_${i}.fits" \
                ifmt=fits \
                omode=meta \
               | tee "datatable_crossmatched_${i}.meta.txt"
               # meta
fi


# 
# Now we cross-match the A3COSMOS catalog
# 
i=$((i+1))
echo "************"
echo "*** ${i} ***"
echo "************"
if [[ ! -f "datatable_crossmatched_${i}.fits" ]]; then
topcat -stilts tmatch2 \
                in1="datatable_crossmatched_$((i-1)).fits" \
                ifmt1=fits \
                values1="RA Dec" \
                suffix1="" \
                join=all1 \
                find=all \
                in2="$input_cat_Liu2020a" \
                icmd2='addcol "WAVELENGTH_ALMA" -units "um" "wObs"' \
                icmd2='addcol "FLUX_ALMA" -units "mJy" "fObs"' \
                icmd2='addcol "FLUXERR_ALMA" -units "mJy" "dfObs"' \
                icmd2='keepcols "ID RA Dec WAVELENGTH_ALMA FLUX_ALMA FLUXERR_ALMA"' \
                values2="RA Dec" \
                suffix2="_Liu2020a" \
                matcher=sky \
                params="${crossmatched_seplimit}" \
                fixcols=dups \
                ocmd='replacecol "GroupID" -name GroupID_Liu2020a "GroupID"' \
                ocmd='replacecol "GroupSize" -name GroupSize_Liu2020a "GroupSize"' \
                ocmd='replacecol "Separation" -name Sep_Liu2020a -units "arcsec" "Separation"' \
                ofmt=fits \
                out="datatable_crossmatched_${i}.fits"
                # 
topcat -stilts tpipe \
                in="datatable_crossmatched_${i}.fits" \
                ifmt=fits \
                omode=meta \
               | tee "datatable_crossmatched_${i}.meta.txt"
               # meta
fi


# 
# Finish
# 
i=$((i+1))
echo "*************"
echo "*** Final ***"
echo "*************"
#if [[ ! -f "datatable_crossmatched_final.fits" ]]; then
echo cp "datatable_crossmatched_$((i-1)).fits" "datatable_crossmatched_final.fits"
cp "datatable_crossmatched_$((i-1)).fits" "datatable_crossmatched_final.fits"
#fi






