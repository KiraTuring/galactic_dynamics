# Galaxy Pipeline Troubleshooting

Common issues across the JAM and Schwarzschild preprocessing pipelines.

## JAM Modeling

### `AssertionError: Inclination too low` / `q_lum < ...`

**Cause:** `q.hi` ≥ MGE `q_min`. The inclination formula requires `q < q_min`.

**Fix:** Set `q.hi` strictly less than `mge.q_min` from `description.yaml`.

### `KeyError: 'lg_ml_inf'`

**Cause:** The JAM PARAMS list requires all 13 parameters, even when
`varying_ml: false`.

**Fix:** Add `lg_ml_inf` to the config with `fixed: true`.

### `AttributeError: 'PosixPath' object has no attribute 'find'`

**Cause:** `use_hist: true` but no histogram file exists.

**Fix:** Set `kinematic_settings.use_hist: false`.

### JAM process dies silently in background

**Cause:** `nohup` alone doesn't fully detach; shell exit kills the child.

**Fix:** Use `setsid` instead of `nohup`:

```bash
setsid python -u scripts/run_jam.py config.yaml > /tmp/jam.log 2>&1 &
```

### Stale `.pyc` after editing `jam_fit/` source

**Cause:** Python bytecode cache retains old function signatures.

**Fix:**
```bash
rm -rf JAM/jam_fit/__pycache__ JAM/__pycache__
```

## PSF Fitting

### `fit_psf.find_hst_file()` raises `FileNotFoundError`

**Cause:** The function searches for `*f555w*_pc_drz.fits` or `*f555w*drz.fits`.
Non-F555W bands (e.g. F702W) won't match.

**Fix:** Either rename the file or monkey-patch:
```python
import fit_psf
fit_psf.find_hst_file = lambda g: f"data/raw/{g}/my_hst_file.fits"
```

### SAURON PSF fit: `tricontour` crashes with "non-finite values"

**Cause:** `read_sauron_kin()` in `jam_fit.prep` computes `-2.5*log10(flux)`
where SAURON C2D flux contains ≤0 values, producing NaN.

**Fix:** Monkey-patch `read_sauron_kin` to sanitize the flux array:

```python
flux_safe = np.where(flux > 0, flux, np.nan)
flux_plot = -2.5 * np.log10(flux_safe / np.nanmax(flux_safe))
valid = np.isfinite(flux_plot)
plt.tricontour(xpix[valid], ypix[valid], flux_plot[valid], ...)
```

### PSF fit hits bounds (sigma1→0.1, sigma2→6.0, chi2 >100)

**Cause:** Poor HST image quality — single-exposure C0F files without
drizzling, or incorrect chip extraction from multi-chip WFPC2 data.

**Fix:** Only use HLA DRZ products (not C0F). If no DRZ available, estimate
PSF from FWHM ÷ 2.355 instead.

### `append_log` crashes with `FileNotFoundError`

**Cause:** `fit_psf.py` tries to append to `data/processed/{galaxy}/description.yaml`
but the file doesn't exist.

**Fix:** Create a stub first:
```bash
echo "galaxy: {galaxy}" > data/processed/{galaxy}/description.yaml
```

## Data Preprocessing

### HST data not found in MAST

**Possible reasons:**
1. Target name has a suffix like `{galaxy}-NUC` — try both
2. Different proposal — check multiple proposal IDs
3. No HST data exists — use FWHM estimate

### SAURON PA mismatch >30°

**Cause:** OASIS and SAURON velocity field PAs differ at large radii.

**Consequence:** Galaxy may NOT be axisymmetric. Exclude from JAM/AxiSchw
axisymmetric models.

### MAST 503 errors

**Cause:** Portal API flaky.

**Fix:** Use direct download URL:
```
https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLA/url/...
```

## Miscellaneous

### Relative paths fail in `run_jam.py`

**Cause:** The config uses `io_settings.input_directory: ../data/...`.
Must run from the `JAM/` directory for relative paths to resolve.

**Fix:** Always `cd JAM` first.
