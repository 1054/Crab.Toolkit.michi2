# NOTES — known code issues found during the 2026-08 documentation pass

Documentation-only branch (`docs-improvement`). These are bugs observed and
verified while testing the demos; **no code was changed on this branch**.
Fixes should go on a separate code branch.

## 1. `michi2-run-SED-fitting-v5` — inverted auto-deploy condition

Around the "Checking libraries:" block:

```bash
if [[ $check_lib_files -lt ${#loop_lib_files[@]} ]] || [[ $Overwrite -ge 3 ]]; then
    michi2-deploy-files ${deploy_params[@]}
else
    echo "All libraries found."
fi
```

`check_lib_files` counts MISSING libraries, so the condition is true only
when **at least one** library already exists. In a completely fresh working
directory (nothing deployed) it is false, no deployment happens, the script
prints the misleading "All libraries found." and `michi2_v05` then fails
with `The input library data "lib.*.SED" does not exist!`.

Suggested fix: invert the test (`-gt 0` instead of `-lt count`), or simply
always deploy (deploy is idempotent apart from re-unzipping).

## 2. `michi2-deploy-files` vs wrapper default stellar library mismatch

`michi2-deploy-files SED` (no extra args) deploys the **BC03.MultiAge**
stellar library (its own default), while `michi2-run-SED-fitting-v5`
defaults to `-lib-stellar BC03.200Myr`. Running the wrapper with defaults
right after a default deploy therefore fails on the missing
`lib...Age200Myr.EBV.SED`. Either align the two defaults or always pass an
explicit `-lib-stellar`.

## 3. `michi2_plot_SED_fitting_results_for_michi2_v05.py` — all-zero
   `best-fit_param_*.txt` under numpy ≥ 2

In `analyze_chisq_distribution()`:

```python
if param_stats['valid'] is True:
```

With numpy ≥ 2, `crab_bin_compute_param_chisq_histogram` returns
`valid` as `numpy.bool_(True)`; the identity test `is True` fails and the
zero branch writes all-zero parameter files. (Under the older numpy used
for the published ID_628146 demo results it was a Python `bool` and the
values were written.) The `chi-square_table_*.txt` files stay valid.

Suggested fix: `if bool(param_stats['valid']):`

## 4. `michi2_filter_flux.py` — IndexError under numpy ≥ 2 on tables with
   duplicated wavelengths

`data_table.remove_rows(...)` is called with row indices collected BEFORE
earlier removals, so a later `data_table[isel2]` slices with stale indices
(`IndexError: index 125 is out of bounds for axis 0 with size 123`). Hit on
the ID_628146 demo datatable (which contains wavelength duplicates).

Suggested fix: recompute indices after each removal (or remove all rows in
one call with a merged, sorted index array).

## 5. `-trial` prints `Trial=Trial`

Cosmetic: `echo Trial=Trial` should be `echo Trial=$Trial`.

---

**Update (2026-08-15):** all five issues above are fixed and merged into
`master` (one commit per bug, each verified against the demos —
empty-directory auto-deploy, default-lib alignment, real values in
`best-fit_param_*.txt` under numpy 2.4, duplicate-wavelength tables
completing).
