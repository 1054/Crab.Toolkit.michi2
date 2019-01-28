#!/bin/bash

# copy file
if [[ ! -f 'datatable_photometry_with_NOEMA_1mm.txt' ]]; then
    cp '../Galsed_z36_sample_with_NOEMA_1mm/datatable_photometry_with_NOEMA_1mm.txt' 'datatable_photometry_with_NOEMA_1mm.txt'
fi


# cross-match with Skelton+2014 GOODS-N 3D-HST photometry catalog
topcat -stilts tmatchn nin=2 \
               in1='datatable_photometry_with_NOEMA_1mm.txt' ifmt1=ascii \
               in2='/Volumes/GoogleDrive/Team Drives/DeepFields/Catalogs/by_Field/GOODSN/goodsn_3dhst.v4.1.cat.FITS' ifmt2=fits \
               matcher=sky \
               params=1 \
               icmd1='addcol ref_z -after z \"zphot_Liu2018\"' \
               values1="ra de" \
               values2="ra dec" \
               fixcols=all \
               suffix1='_Liu2018' \
               suffix2='_3DHST' \
               out='datatable_photometry_with_NOEMA_1mm_with_optical.fits'

echo "Output to 'datatable_photometry_with_NOEMA_1mm_with_optical.fits'"




