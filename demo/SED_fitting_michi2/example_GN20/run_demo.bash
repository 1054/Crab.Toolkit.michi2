#!/bin/bash
#
# Minimal self-contained michi2 SED-fitting example — GN20 (z = 4.055).
#
# Photometry credit: GN20 (a FIR-luminous submillimeter galaxy at z=4.055),
# 10-band subset (optical i-band to SPIRE 500um) converted to mJy from the
# SED_Fitting_Expert GN20 fixture, which itself derives from the observed
# flux table published at https://zenodo.org/record/3958272 (Liu+2021).
#
# This script completes in a few minutes on a laptop:
#   - step 1 deploys the (zipped) model libraries once;
#   - step 2 is a quick smoke test (-trial, seconds);
#   - step 3 is a science-grade run (-sampling 3000, a few minutes).
#
# Usage (BASH shell only):
#     source /path/to/Crab.Toolkit.michi2/SETUP.bash
#     bash run_demo.bash
#

set -e

FLUX_FILE="GN20_flux.txt"
REDSHIFT="4.055"

# sanity check
if [[ ! -f "$FLUX_FILE" ]]; then
    echo "Error! \"$FLUX_FILE\" not found in $(pwd)!"
    exit 255
fi

# 1. Deploy the model libraries into this directory (required once;
#    the wrapper's auto-deploy does not trigger in a completely empty dir).
michi2-deploy-files SED

# 2. Quick smoke test (sampling 30 — checks the pipeline only; the fitted
#    parameter files from a trial run are NOT meaningful).
michi2-run-SED-fitting-v5 -flux "$FLUX_FILE" -redshift $REDSHIFT \
    -lib-stellar BC03.MultiAge -obj-name GN20 -trial -overwrite

# 3. Science-grade run (a few minutes with 4 cores; use -sampling 150000
#    for publication-level runs — takes hours).
michi2-run-SED-fitting-v5 -flux "$FLUX_FILE" -redshift $REDSHIFT \
    -parallel 4 -sampling 3000 \
    -lib-stellar BC03.MultiAge -lib-dust DL07UPD2010 \
    -obj-name GN20 -overwrite

echo ""
echo "Done. Key outputs:"
echo "  fit_5.out                              (chi2 + params per combination)"
echo "  results_fit_5/best-fit_param_Mstar.txt (log10 Msun statistics)"
echo "  results_fit_5/best-fit_SED_GN20.txt    (best-fit SED curve, obs-frame um + mJy)"
echo "  results_fit_5/fit_5.pdf                (best-fit SED figure)"
