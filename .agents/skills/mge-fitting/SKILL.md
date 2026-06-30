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
conda activate schw   # mgefit + astropy + matplotlib
```

Required data:

| Data | Location | Format |
|------|----------|--------|
| HST WFPC2 image | `data/raw/{galaxy}/hst_*f814w*_pc_drz.fits` | HLA drizzled FITS (SCI extension) |
| Same, F555W fallback | `data/raw/{galaxy}/hst_*f555w*_pc_drz.fits` | Same structure |

Filter is auto-detected from `FILTNAME1` header keyword. Supported filters:
F555W (4.82), F606W (4.77), F702W (4.71), F814W (4.56), F160W (5.36).

## Pipeline

```
HST drz FITS  ──►  sky subtraction (auto/user)
              ──►  find_galaxy()          ──► center, eps, PA
              ──►  sectors_photometry()   ──► radial SB in sectors
              ──►  mge_fit_sectors_regularized()  ──► total_counts, sigma, q
              ──►  unit conversion        ──► L⊙ pc⁻²
              ──►  save ECSV + diagnostics
```

## Script

`scripts/fit_mge_hst.py`

```bash
python scripts/fit_mge_hst.py <galaxy> [options]
```

### Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `galaxy` | (required) | Galaxy name (e.g. NGC4621, NGC5845) |
| `--fwhm` | `0.13` | PSF FWHM in arcsec (WFPC2 F814W default; F555W ~0.10) |
| `--ngauss` | `20` | Initial number of Gaussian components |
| `--trim` | `0.05` | Fractional edge trim |
| `--sky` | auto | Sky background in counts/s/pix (auto-estimate from corners) |

### Output (`data/processed/{galaxy}/`)

File names use the actual HST filter (e.g. `mge_F814W.ecsv`, `mge_F555W.ecsv`).

| File | Content |
|------|---------|
| `mge_{filter}.ecsv` | MGE parameters (I, sigma, q, pa_twist) |
| `mge_{filter}_components.png` | Gaussian ellipses + radial profile |
| `mge_{filter}_contours.png` | Two-panel contour: full FOV + zoom 3″×3″ (black=data, red=model) |
| `mge_{filter}_radial.png` | Radial profile (major + minor axis with data points) |
| `mge_{filter}_comparison.png` | Comparison with existing `mge.ecsv` (if found) |

### Unit conversion

```
counts/s/pix  →  PHOTFLAM  →  f_λ (erg/cm²/s/Å)
                           →  f_ν = f_λ × λ² / c
                           →  μ_AB = -2.5 log₁₀(f_ν) - 48.6
                           →  I (L⊙ pc⁻²) = (206265)² × 10^(-0.4 × (μ_AB - M⊙))
```

Solar AB magnitude `M⊙` is looked up by filter (Willmer 2018):
- F555W: 4.82, F814W: 4.56

## Coordinate convention

`mgefit` internally uses `(row, col)` convention matching numpy array indexing
(`img[xc, yc]`).  The same convention is used throughout the pipeline —
**no coordinate swapping is needed**.

```python
sec = find_galaxy(img, ...)                         # sec.xmed = row, sec.ymed = col
pho = sectors_photometry(img, eps, theta, sec.xmed, sec.ymed, ...)
```

## Number of components

The regularized fit (`mge_fit_sectors_regularized`, Cappellari 2002) eliminates
Gaussians that are not statistically required.  Typical counts for WFPC2 PC
images (~24″ radius field):

| Galaxy | Components | Note |
|--------|-----------|------|
| NGC4621 (Reff ~35″) | 5 | Galaxy larger than field |
| NGC5845 (compact) | 7 | Fits better in PC field |

For a full MGE covering all radii, combine PC + WF + ground-based data
(Cappellari et al. 2006, MNRAS, 366, 1126 — the SAURON project).

## Common issues

### Too few Gaussian components

The regularized fit drops unneeded Gaussians.  To get more:
- Increase `--ngauss` (e.g. `--ngauss 30`)
- Combine with wider-field data to better constrain the outer profile

### Surface brightness mismatch with existing MGE

Existing MGE catalogs (SAURON/ATLAS3D) use ground-based SDSS data with
different filters, seeing, and FOV.  A factor of 10–100× in central surface
brightness between HST and ground-based MGEs is expected and normal.

### Contours don't match (offset or rotation)

- Ensure `sigmapsf` was passed correctly (the `--fwhm` parameter)
- The script now plots contours manually (avoids `mge_print_contours` bugs)

### Negative / bad pixels in image

The script reads WHT extension and masks pixels with `WHT <= 0`.  Increase
`--trim` if edge artifacts remain.

## Scripts reference

| Script | Location | Purpose |
|--------|----------|---------|
| `fit_mge_hst.py` | `.agents/skills/mge-fitting/scripts/` | MGE fitting from HST images |

## References

- Cappellari (2002, MNRAS, 333, 400) — MGE fitting method
- Cappellari et al. (2006, MNRAS, 366, 1126) — SAURON project MGE models
- Willmer (2018, ApJS, 236, 47) — Solar absolute AB magnitudes
