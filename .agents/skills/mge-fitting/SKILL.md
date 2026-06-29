---
name: mge-fitting
description: >
  Fit Multi-Gaussian Expansion (MGE) photometric models from HST images
  using mgefit. Covers HST image selection, sector photometry, MGE fitting
  with regularization, unit conversion (counts → L⊙ pc⁻²), and diagnostic
  plotting. Use whenever creating MGE models for JAM or Schwarzschild
  dynamical modeling from HST imaging data, or when replacing an existing
  MGE with a higher-resolution HST-based fit.
---

# MGE Fitting from HST Images

Fit a Multi-Gaussian Expansion (MGE) photometric model to an HST WFPC2 image
using the `mgefit` package (Cappellari 2002). The MGE decomposes the galaxy's
surface brightness into a sum of 2D Gaussians, which is the required input
for JAM and Schwarzschild dynamical models.

## When to use

- Creating an MGE model from an existing HST image
- Replacing an existing MGE with a higher-resolution HST-based fit
- Adding a new galaxy to the JAM or AxiSchw modeling pipeline
- Testing PSF and MGE systematics with different HST filters

## Prerequisites

```bash
conda activate schw   # mgefit is installed in this env
```

Required data:

| Data | Location | Format |
|------|----------|--------|
| HST WFPC2 F814W image | `data/raw/{galaxy}/hst_*f814w*_pc_drz.fits` | HLA drizzled FITS (SCI extension) |
| HST WFPC2 F555W image | `data/raw/{galaxy}/hst_*f555w*_pc_drz.fits` | Same format (alternative filter) |

The script reads the `SCI` extension and uses `PHOTFLAM` / `CENTRWV` header
keywords for photometric calibration.

## Pipeline

```
HST F814W drz FITS  ──►  find_galaxy()      ──► center, eps, PA
                             │
                             ▼
                     sectors_photometry()   ──► radial SB in sectors
                             │
                             ▼
                     mge_fit_sectors()      ──► total_counts, sigma, q
                             │
                             ▼
                     Unit conversion        ──► L⊙ pc⁻²
                             │
                             ▼
                     Save + diagnostics     ──► mge_f814w.ecsv, plots
```

## Script

`scripts/fit_mge_hst.py`

```bash
python scripts/fit_mge_hst.py <galaxy> [--fwhm FWHM] [--ngauss N]
```

### Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `galaxy` | (required) | Galaxy name (e.g. NGC4621) |
| `--fwhm` | `0.13` | PSF FWHM in arcsec (WFPC2 F814W default) |
| `--ngauss` | `15` | Initial number of Gaussian components |
| `--trim` | `0.05` | Fractional edge trim |

### Output (`data/processed/{galaxy}/`)

| File | Content |
|------|---------|
| `mge_f814w.ecsv` | MGE parameters (I, sigma, q, pa_twist) |
| `mge_f814w_components.png` | Gaussian components visualization |
| `mge_f814w_contours.png` | Data (black) vs model (red) contour overlay |
| `mge_f814w_radial.png` | Radial profile fit |
| `mge_f814w_comparison.png` | Comparison with existing MGE (if `mge.ecsv` exists) |

### Unit conversion

```
counts/s/pixel  →  PHOTFLAM  →  f_λ (erg/cm²/s/Å)
                             →  f_ν = f_λ × λ² / c
                             →  μ_AB = -2.5 log₁₀(f_ν) - 48.6
                             →  I (L⊙ pc⁻²) = (206265)² × 10^(-0.4 × (μ_AB - M⊙_F814W))
```

Solar AB magnitude in F814W: M⊙ = 4.56 (Willmer 2018).

## MGE coordinate convention

`mgefit` uses a non-standard coordinate convention internally. Both
`sectors_photometry` and `mge_print_contours` — as well as `find_galaxy` —
treat the first spatial argument (`xc`) as the **row** index and the second
(`yc`) as the **column** index (i.e., `img[xc, yc]` in numpy). This differs
from the astronomical convention where X = column, Y = row.

**No swapping is needed** — just pass `sec.xmed, sec.ymed` directly:

```python
# Correct — x=row, y=col throughout
sec = find_galaxy(img, ...)
pho = sectors_photometry(img, sec.eps, sec.theta, sec.xmed, sec.ymed, ...)
mge_print_contours(img, sec.theta, sec.xmed, sec.ymed, sol, ...)
```

## Number of components

The regularized fit (`mge_fit_sectors_regularized`) eliminates Gaussians that
are not statistically required. For HST PC images (~24″ radius field), the
number of surviving components depends on how much of the galaxy profile is
covered:

| FOV | Typical components | Note |
|-----|-------------------|------|
| WFPC2 PC (~24″) | 2–5 | Galaxy fills most of field |
| WFPC2 WF (~150″) | 5–8 | Better outer constraints |
| Ground-based (~5′) | 10–17 | Full galaxy profile |

If the galaxy's effective radius exceeds the image field, the MGE model is
only reliable for the inner region. Combine HST and ground-based data for a
full MGE covering all radii.

## Common issues

### Contour center doesn't match data center

Check coordinate conventions (see above). Ensure `sigmapsf` and `normpsf` are
passed to `mge_print_contours`:

```python
mge_print_contours(img, theta, xc, yc, sol,
                   sigmapsf=sigma_psf, normpsf=[1.0], scale=pixscale)
```

### Too few Gaussian components

The regularized fit drops unneeded Gaussians. If more components are desired:

- Increase `--ngauss` (e.g. `--ngauss 30`)
- Remove `outer_slope` constraint (currently always applied by mgefit)
- Combine with wider-field data to better constrain the outer profile

### Surface brightness mismatch with existing MGE

Existing MGE catalogs (e.g. SAURON/ATLAS3D) are typically based on ground-based
SDSS data with different filters, seeing, and photometric calibration. A factor
of 10–100× difference in central surface brightness is expected between HST and
ground-based MGEs.

### Negative pixel values in image

Trim edge pixels or apply a brighter `--trim` threshold. The sector photometry
ignores pixels below `minlevel`.

## Scripts reference

| Script | Location | Purpose |
|--------|----------|---------|
| `fit_mge_hst.py` | `.agents/skills/mge-fitting/scripts/` | MGE fitting from HST images |

## References

- Cappellari (2002, MNRAS, 333, 400) — MGE method
- Scott et al. (2013, MNRAS, 432, 1894) — ATLAS3D MGE catalog
- Willmer (2018, ApJS, 236, 47) — Solar absolute magnitudes
