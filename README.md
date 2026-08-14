# Crab.Toolkit.michi2

**michi2** is a _Crab_ toolkit for minimizing chi-square in any 1D model fit —
for example, galaxy spectral energy distribution (SED) fitting and molecular
line ladder (SLED) fitting. The SED fitting models a galaxy as a linear
combination of up to five components — BC03 stellar populations, Mullaney+2011
AGN, Draine & Li (2007) warm and cold dust, and a radio power law — and
exhaustively loops over all model-library combinations, giving fine grids of
ISRF (Umin) and dust mass.

> **If you use this code in a publication, please cite**
> [Liu 2020, ascl:2005.002](https://ui.adsabs.harvard.edu/abs/2020ascl.soft05002L)
> (the code paper) and
> [Liu et al. 2021, ApJ, 909, 56](https://ui.adsabs.harvard.edu/abs/2021ApJ...909...56L)
> (the science application). BibTeX at the bottom of this page.

**Requirements:** bash (macOS or Linux) · a Python 3 with
`numpy`/`matplotlib`/`astropy` on the `PATH` (for post-processing and
plotting) · pre-compiled `michi2_v05` binaries ship in `bin/` (macOS x86_64,
Linux glibc 2.14 / 2.22). The C source lives in a separate repository:
<https://gitlab.com/1054/michi2>.

---

## Quickstart (SED fitting)

```bash
# 1. download/clone this repo, then (in a BASH shell):
source /path/to/Crab.Toolkit.michi2/SETUP.bash

# 2. in an EMPTY working directory, deploy the model libraries once:
michi2-deploy-files SED

# 3. prepare a 3-column flux file (see "Input format" below), then:
michi2-run-SED-fitting-v5 -redshift 1.5 -flux extracted_flux.txt -parallel 2

# outputs: fit_5.out and results_fit_5/ (parameters, best-fit SED, PDFs)
```

A ready-made runnable example is in
[`demo/SED_fitting_michi2/example_GN20/`](demo/SED_fitting_michi2/example_GN20/)
(completes in minutes). A full A³COSMOS demo is in
[`demo/SED_fitting_michi2/ID_628146/`](demo/SED_fitting_michi2/ID_628146/).

---

# Usage for SED Fitting #

The SED fitting contains 5 components (or any number of components): BC03
stellar, AGN, DL07 warm dust, DL07 cold dust, and radio. The radio component
is constrained by the DL07 warm+cold summed rest-frame 8–1000 µm luminosity
through the IR–radio correlation (`-qIR`, default q = 2.4), so it is not
independent. If you do not want that, simply do not include radio data in
the input file (or pass `-free-radio`).

The fitting needs a specific redshift as the input, so it **cannot fit the
redshift**.

The advantage of this code is that it loops over all DL07 dust models and
combinations with AGN models and stellar models, so it can fit ISRF and
Mdust in a finer parameter grid than template-averaging codes. The downside
is speed: a full-sampling 5-component fit takes minutes to hours (see
"-sampling" below; use `-parallel N` on multi-core machines).

## 1. Download the code

Scroll up, there is a green button **"Clone or download"**. Click it and
select **"Download ZIP"**, then uncompress it somewhere, e.g. your working
directory. Currently the ZIP file size is about 52 MB. If you are familiar
with git you can also clone it, but note the repository size is a bit large
(hundreds of MB).

## 2. Installation

No compilation is needed for normal use — pre-compiled binaries ship with
the repo:

| Platform | Binary | Note |
|----------|--------|------|
| macOS | `bin/bin_mac/michi2_v05_mac` | x86_64; runs on Apple Silicon via Rosetta 2 (verified) |
| Linux (glibc ≥ 2.14) | `bin/bin_linux_glibc_2_14/michi2_v05_linux_glibc_2_14` | x86_64 |
| Linux (glibc ≥ 2.22) | `bin/bin_linux_glibc_2_22/…` | x86_64 |

`bin/michi2_v05` is a small dispatch script that picks the right binary per
OS. The **C source code lives in a separate repository**:
<https://gitlab.com/1054/michi2> — clone it there if you need to rebuild or
modify the core fitter.

Setup (must be under a BASH shell, not zsh/fish):

```bash
bash
source /some/path/Crab.Toolkit.michi2/SETUP.bash
```

This prepends `<repo>/bin` to your `PATH`. The wrapper scripts also
self-source `bin/bin_setup.bash` when needed, so an absolute-path call works
too.

The post-processing/plotting step calls Python scripts that import
`numpy`, `matplotlib` and `astropy` — make sure a Python 3 with those
packages is on your `PATH`.

## 3. Deploy the model libraries (REQUIRED once per working directory)

The model libraries ship **zipped** in `data/lib_SED/*.zip` and must be
extracted into your working directory before the first run, together with
the `filters/` directory and the `filter.list` file (61 filter curves):

```bash
cd /path/to/your/working/directory/
michi2-deploy-files SED            # deploy the defaults:
                                   #   BC03.MultiAge + MullaneyAGN +
                                   #   DL07UPD2010 + Radio + filters
michi2-deploy-files SED BC03.400Myr   # ...or a specific stellar library
```

`michi2-run-SED-fitting-v5` nominally deploys missing libraries itself, but
its auto-deploy check currently has an inverted condition: in a **completely
empty** directory it prints "All libraries found." and the fit then fails
with `The input library data "lib.*.SED" does not exist!`. **Always run
`michi2-deploy-files` yourself in a fresh directory** (the wrapper's
auto-deploy only helps when at least one library file is already present).

⚠️ **Library-selection caveat:** `michi2-deploy-files SED` (no extra args)
deploys the **BC03.MultiAge** stellar library, but the wrapper's default
`-lib-stellar` is **BC03.200Myr**. Either pass the same explicit
`-lib-stellar` to the wrapper as what you deployed, or deploy matching
libraries, e.g.:

```bash
michi2-deploy-files SED BC03.200Myr    # then run WITHOUT -lib-stellar
# or
michi2-deploy-files SED                # then run with -lib-stellar BC03.MultiAge
```

Re-running `michi2-deploy-files` re-extracts the libraries (it moves the
previous `libs.list` to `libs.list.backup`); `filter.list` is never
overwritten once it exists.

## 4. Input format

Each source needs one SED fitting input text file with **three columns**
separated by white space:

1. wavelength in micron (observed frame),
2. flux density in milli-Jansky,
3. flux-density error in milli-Jansky.

Rows can be commented out with `#`. Use S/N > 3 bands (upper limits are not
supported by the default fitting path — drop them). Two extra columns (flux
unit, filter name) are ignored by the fitter but are handy for bookkeeping
and are echoed by the plotting script.

Example (`extracted_flux.txt`):

```
# wavelength_um     flux_mJy    flux_err_mJy
    3.56343       0.12417548     0.012417548
    4.51101       0.09488227     0.009488227
    5.75934      0.075890755    0.0075890755
    7.95949       0.08685837     0.009947749
       24.0        1.3679016      0.13679016
      100.0        41.577702       4.1577702
      160.0          60.0042       6.3393602
      250.0          58.0821         5.80821
      350.0        36.100101       6.0408401
      500.0          18.878901       2.2385001
```

Note the wavelengths should be unique per file (duplicated wavelengths can
confuse the preparation scripts).

## 5. Run michi2

Having `source`d `SETUP.bash` and deployed the libraries, change to the
directory with your flux file and run:

```bash
michi2-run-SED-fitting-v5          # without arguments prints the usage

michi2-run-SED-fitting-v5 -redshift 1.5 -flux "extracted_flux.txt" -parallel 2

# a full-featured invocation:
michi2-run-SED-fitting-v5 -redshift 1.5 -flux "extracted_flux.txt" \
    -parallel 2 -sampling 150000 \
    -lib-stellar BC03.MultiAge -lib-dust DL07UPD2010 \
    -lib-AGN MullaneyAGN -lib-radio Radio \
    -freeze-radio -qIR 2.4 -Umin 1 -minEBV 0.2 \
    -obj-name "My Galaxy" -overwrite
```

The SED fitting is **slow at full sampling** — a 5-component fit easily
takes hours on a laptop because it loops over all combinations of all input
models. It is parallelized (`-parallel N`); for a first look use `-trial`
or `-sampling 3000`.

### Command-line reference (`michi2-run-SED-fitting-v5`)

**Required:**

| Flag | Meaning |
|------|---------|
| `-flux FILE` | the 3-column input file (§4) |
| `-redshift Z` | the fixed redshift (**cannot be fitted**) |

**Runtime:**

| Flag | Default | Meaning |
|------|---------|---------|
| `-parallel N` | 2 | number of CPU cores |
| `-sampling N` | 10000 | parameter-space sampling density. 30 ≈ instant smoke test (`-trial` sets this), 3000 ≈ science-grade for a ~20-band SED, 150000 ≈ the value used in Liu et al. (2021) |
| `-trial` | — | fast trial run (forces sampling 30, parallel 2) |
| `-overwrite` | off | re-fit even if `fit_5.out` exists. It is a **counter**: pass once to re-fit; twice to also rebuild `fit_5.in` from `-flux` |
| `-fit-name NAME` | `fit_5` | output prefix (`NAME.out`, `results_NAME/`) |
| `-obj-name NAME` | `obj` | source name; names `best-fit_SED_<NAME>.txt` and labels plots |
| `-filter-flux` | off | run the flux-filtering scripts (drop low-S/N bands) when building `fit_5.in` |

**Library selection:**

| Flag | Accepted values |
|------|-----------------|
| `-lib-stellar` | `BC03`, `BC03.200Myr` (wrapper default), `BC03.400Myr`, `BC03.MultiAge`, `FSPS.CSP.tau.0p1Gyr`, `FSPS.CSP.tau.1Gyr` |
| `-lib-dust` | `DL07`, `DL07UPD2010` (default), `DL07UPD2010FIR40122` |
| `-lib-AGN` | `MullaneyAGN` |
| `-lib-radio` | `Radio` |

**Component toggles** (default: all five components on):

| Flag | Meaning |
|------|---------|
| `-no-stellar` / `-no-AGN` / `-no-dust` / `-no-radio` | drop that component. `-no-dust` automatically sets `-no-radio` and disables `-freeze-radio` |
| `-no-warm-dust` / `-no-cold-dust` | use only one of the two DL07 components (cannot set both) |

**Physics constraints:**

| Flag | Default | Meaning |
|------|---------|---------|
| `-freeze-radio` | off | output a radio component frozen to the IR–radio relation (with `-qIR`) even when no radio data are fit |
| `-free-radio` | off | do NOT tie the radio normalization to the dust luminosity |
| `-free-dust` | off | decouple the warm and cold dust library indices (by default they share one DL07 grid point / Umin) |
| `-qIR Q` | 2.4 | the q IR–radio correlation used to tie radio to L_IR(8–1000 µm) |
| `-Umin X` | — | lower bound on the dust Umin grid |
| `-qPAH X` | — | upper bound on qPAH |
| `-minEBV X` | — | lower bound on the stellar E(B−V) (BC03 libraries only) |

## 6. Output reference

Everything is written to the current working directory:

```
fit_5.out          # one row per library-combination trial: indices,
                   # normalizations (a_i), chi2, per-library physical params
fit_5.out.info     # basic info for later steps
results_fit_5/     # everything below
```

### `results_fit_5/best-fit_param_<NAME>.txt` — fitted parameters

One file per parameter, with five statistics — `param_median`,
`param_best`, `param_sigma`, `param_L68`, `param_H68` — where `best` is the
value at the minimum chi-square and median/L68/H68 come from the Δχ² = 5.89
window (5-parameter 1σ). The file has three header lines and a `null`
template row before the data row (take the **last** row with five numbers).

Units — **pay attention to which parameters are log10**:

| Parameter | Unit | Scale |
|-----------|------|-------|
| `Mstar` | M☉ | **log10** |
| `Age` | Gyr | **log10** |
| `EBV` | mag | linear |
| `LAGN` | L☉ | **log10** |
| `LIR_warm`, `LIR_cold`, `LIR_total` (8–1000 µm), `LIR_total_40_400` (40–400 µm) | L☉ | **log10** |
| `Mdust_warm`, `Mdust_cold`, `Mdust_total` | M☉ | **log10** |
| `Umin_warm`, `Umin_cold`, `Umean_total` | — | linear (DL07 radiation field) |
| `fPDR_total` | fraction | **log10** (see "Definitions & caveats") |
| `qPAH_warm`, `qPAH_cold` | — | linear |
| `AGN_TYPE` | — | Mullaney template index (linear) |

`results_fit_5/chi-square_table_<NAME>.txt` holds the raw (chi2, param)
pairs behind these statistics (params stored **linear** there; michi2
log-transforms when computing the stats).

### `results_fit_5/best-fit_SED_<obj-name>.txt` — best-fit SED curve

Two columns: `wavelength_obsframe_um`, `flux_obsframe_mJy` — 40001
log-spaced points in the **observed frame** (already redshifted; do not
multiply by 1+z). This is the total (all components summed); per-component
curves can be rebuilt from the `fit_5.out` normalizations plus the library
files (that is what `fit_5.pdf` plots).

### Plots

| File | Content |
|------|---------|
| `results_fit_5/fit_5.pdf` | best-fit SEDs with component decomposition |
| `results_fit_5/fit_5.best.pdf` | the single best-fit SED |
| `results_fit_5/fit_5.chisq.pdf`, `fit_5.best.chisq.pdf` | chi-square histograms per parameter |

## 7. Optionally re-plotting chi2 distribution and compute best-fits

Assuming you have already `source`d the `SETUP.bash`:

```bash
michi2-plot-fitting-results          # without arguments prints the usage

michi2-plot-fitting-results fit_5.out -flux extracted_flux.txt -source YOUR_SOURCE_NAME
```

## 8. Definitions & caveats

- **fPDR is a mass fraction, not DL07 Eq. 29's luminosity fraction.** Our
  fPDR is the mass fraction of Umin–Umax dust — i.e. it is **DL07 γ**
  (M_dust(PDR) / M_dust(total)) — *not* the f_PDR of Draine & Li (2007)
  Eq. 29, which is the luminosity fraction of U > 100 dust. Keep this in
  mind when comparing with DL07-based codes.
- **Fixed redshift only.** The redshift is a required input; michi2 does
  not sample it.
- **Upper limits are not fitted** by the default path — provide S/N > 3
  detections only.
- **Known issue (numpy ≥ 2):** with recent numpy versions,
  `best-fit_param_*.txt` files may be written as all zeros. The chi-square
  analysis returns `valid` as `numpy.bool_`, and
  `michi2_plot_SED_fitting_results_for_michi2_v05.py` tests it with an
  identity check (`is True`), which fails → the zero branch is taken. The
  `chi-square_table_*.txt` files remain valid. A one-line fix is to use
  `if bool(param_stats['valid']):`.
- **Known issue (numpy ≥ 2):** `michi2_filter_flux.py` can crash with an
  `IndexError` on tables with duplicated wavelengths (`remove_rows` followed
  by stale-index slicing). Deduplicate wavelengths in the input, or build
  `fit_5.in` yourself (the wrapper skips that step when `fit_5.in` already
  exists and `-overwrite` was passed fewer than two times).

---

# Usage for LVG Fitting #

## LVG Fitting ##

### 1. Get the code (with git)

Same as for SED fitting — download or clone this repository (see §1 above).

### 2. Source the code (must under BASH shell)

```bash
bash
source /some/path/Crab.Toolkit.michi2/SETUP.bash
```

### 3. Prepare line flux data

An example line flux data table, assuming it is named `flux_co_ci.txt`:

```
# X_species  S_species  E_S_species   Molecule
#
101001000  0.21       0.05          CO
101002001  1.0        0.25          CO
101004003  1.68       0.1           CO
101005004  2.2        0.7           CO
101006005  1.8        0.2           CO
101007006  2.19       0.184         CO
102001000  0.70       0.11          C_atom
102002001  1.72       0.20          C_atom
```

The first column is a number needed by our fitting, which is unique for
each line. The first three digits, 101 means CO, and 102 means C_atom. The
second three digits mean the upper level, and the third three digits the
lower level. For example CO J=1-0 is 101 001 000, and CO J=9-8 is 101 009
008. For [CI] it's the same: [CI] 3P1-3P0 is 102 001 000, and [CI] 3P2-3P1
is 102 002 001. The second column is the integrated line flux in Jy km/s,
the third its error. The fourth column is optional; the fitting code only
reads the first three columns.

### 4. Prepare molecular gas Large-Velocity-Gradient model

We also need a molecular gas Large-Velocity-Gradient (LVG) model file before
our fitting. Because the Cosmic Microwave Background (CMB) temperature is
different at different redshift, such a model file needs to be generated for
each redshift.

Try to find a corresponding LVG model file under `data/lib_LVG/`, usually
named like `lib_z_1.500_with_CO_and_C_atom_dV_50.lvg`; copy it to your
working directory and uncompress it if it is a `*.zip` file. Please contact
us for a LVG model file.

(Update on 2022-12-21: for the
[Liu+2022b NGC1365](https://arxiv.org/abs/2212.09652) work I used
`lib_z_0.000_with_CO_and_C_atom_dV_45_vary_XCICO_nonuniform.lvg.zip`, which
allows fitting the [CI/CO] abundance ratio on a non-uniform grid for the
NGC1365 nuclear starburst ring.)

### 5. Run michi2

Having `source`d the `SETUP.bash`, and with the line flux file and LVG model
file in the working directory:

```bash
ls "flux_co_ci.txt"                                          # line flux file
ls "lib_z_4.055_with_CO_and_C_atom_dV_50.lvg"                # LVG model file

michi2-run-fitting-5-components-applying-evolving-qIR        # usage

# one-component fit with 2 CPU cores, sampling 15000 chi-squares
michi2_v05 -obs "flux_co_ci.txt" \
           -lib "lib_z_4.055_with_CO_and_C_atom_dV_50.lvg" \
           -out "result/fit.out" \
           -sampling 15000 \
           -parallel 2 \
           | tee log.txt

# result plots
michi2_plot_LVG_fitting_results.py "result/fit.out"
ls "result/fit.pdf" "result/fit.chisq.pdf" "best-fit_param_"*.txt

# two-component fit: input "-lib" twice and constrain one component
# to always have a lower temperature.
michi2_v05 -obs "flux_co_ci.dat" \
           -lib "lib_z_4.055_with_CO_and_C_atom_dV_50.lvg" \
           -lib "lib_z_4.055_with_CO_and_C_atom_dV_50.lvg" \
           -out "result_two_component_fit/fit.out" \
           -sampling 15000 \
           -parallel 2 \
           -constraint "LIB1_PAR2 < LIB2_PAR2" \
           | tee log_two_component_fit.txt

# result plots
michi2_plot_LVG_fitting_results.py "result_two_component_fit/fit.out"
```

If you have enough lines, say more than 4 lines, it is better to use a
two-component fit.

### 6. Make a nicer CO CI SLED figure

Copy all the files under `demo/LVG_fitting_michi2/plot_nicer_CO_SLED` except
the `result_two_component_fit` subfolder, `flux_co_ci.txt` and
`Plot_nicer_SLED.pdf` to your working directory, then modify
`a_dzliu_code_plot_nicer_SLED.py` to read your own fitting result folder and
line flux data. The output will be your own version of `Plot_nicer_SLED.pdf`.

---

# Citation

If you use michi2, please cite:

```bibtex
@MISC{2020ascl.soft05002L,
  author = {{Liu}, Daizhong},
  title = "{michi2: SED and SLED fitting tool}",
  year = 2020,
  journal = {Astrophysics Source Code Library, record ascl:2005.002},
  url = {https://ascl.net/2005.002}
}

@ARTICLE{2021ApJ...909...56L,
  author = {{Liu}, Daizhong and {Daddi}, Emanuele and {Schinnerer}, Eva and others},
  title = "{CO Excitation, Molecular Gas Density, and Interstellar Radiation
           Field in 76 Local Star-forming Galaxies Detected in CO(2--1) and CO(5--4)}",
  journal = {The Astrophysical Journal},
  year = 2021,
  volume = 909,
  number = 1,
  pages = {56},
  doi = {10.3847/1538-4357/abd801}
}
```
