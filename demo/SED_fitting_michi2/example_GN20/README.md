# Example 1 — minimal michi2 SED fit: GN20 (z = 4.055)

A self-contained, minutes-scale example: 10-band optical-to-FIR photometry
of GN20 (a FIR-luminous submillimeter galaxy at z = 4.055), fitted with all
five components (BC03.MultiAge stellar + Mullaney AGN + DL07 warm/cold dust
+ radio frozen to the IR–radio relation, qIR = 2.4).

Photometry credit: GN20 observed flux table published at
<https://zenodo.org/record/3958272> (Liu et al. 2021), converted to the
michi2 input format (µm / mJy / mJy) via the SED_Fitting_Expert GN20
fixture.

## Run it

```bash
source /path/to/Crab.Toolkit.michi2/SETUP.bash     # BASH shell only
bash run_demo.bash
```

The script (1) deploys the zipped model libraries into this directory,
(2) runs a seconds-scale `-trial` smoke test, and (3) runs a
`-sampling 3000` fit (a few minutes with 4 cores). Note the trial run's
parameter files are not meaningful (sampling 30 only tests the pipeline),
and 3000 is a quick-look density — the library-combination sampling is
stochastic, so use `-sampling 15000+` for reported values (150000 for
publication-level runs).

## Expected outputs

```
fit_5.out                              # chi2 + params per model combination
results_fit_5/best-fit_param_*.txt     # per-parameter statistics
results_fit_5/chi-square_table_*.txt   # raw (chi2, param) pairs
results_fit_5/best-fit_SED_GN20.txt    # best-fit SED curve (obs-frame um + mJy)
results_fit_5/fit_5.pdf                # best-fit SED + component decomposition
```

Reference values from a verified run (`-sampling 3000`, 4 cores; because
the library-combination sampling is stochastic, expect some run-to-run
scatter — the values below are indicative, and use `-sampling 150000` for
publication-level results):

| Parameter | Value | Unit |
|-----------|-------|------|
| Mstar | 11.22 | log10 M☉ |
| LIR_total (8–1000 µm) | 13.18 | log10 L☉ |
| Mdust_total | 9.40 | log10 M☉ |
| EBV | 0.3 | mag |
| Umin (warm = cold) | 40 | — |
| chi2_min | 25.3 | (10 bands) |

Compare with `results_fit_5/fit_5.pdf`: the model passes through the
SPIRE points and shows the warm/cold dust + stellar decomposition.

> **numpy ≥ 2 note:** on unpatched checkouts the `best-fit_param_*.txt`
> files are written as all zeros (a known issue — see NOTES.md #3 / the
> README's "Known issues"; fixed on the `fix-bugs` branch). The numbers
> above were computed from the `chi-square_table_*.txt` files, which stay
> valid either way.
