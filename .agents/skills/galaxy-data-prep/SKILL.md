---
name: galaxy-data-prep
description: >
  Full data preprocessing pipeline for galactic dynamics modeling. Use whenever
  adding a galaxy to the JAM or AxiSchw Schwarzschild pipeline, preprocessing raw
  OASIS/ATLAS3D IFU kinematics, fitting MGE photometric models, downloading HST
  images from MAST, running multi-Gaussian PSF fitting, creating JAM config YAMLs,
  or setting up a new galaxy directory under data/processed/. Covers the full
  workflow from raw FITS to running JAM models.
---

# Galaxy Data Preparation Pipeline

Converts raw OASIS/ATLAS3D integral-field kinematics into JAM-modelling-ready
inputs, with optional HST image download and PSF fitting.

## When to use

- Adding a **new galaxy** to the JAM or AxiSchw modelling pipeline
- **Batch processing** multiple OASIS galaxies from the McDermid+2006 sample
- **Re-fitting the PSF** after downloading new HST images
- **Re-running JAM** with different config parameters after updating PSF/MGE

## Prerequisites

Activate the `schw` conda environment before running any script:

```bash
source /home/haitong/miniconda3/bin/activate schw
```

Required data per galaxy:

| Data | Location | Required? |
|------|----------|-----------|
| OASIS kinematics | `data/raw/{galaxy}/oasis/kinematics_oasis_{galaxy}.fits.gz` | Yes |
| MGE decomposition | `data/processed/sauron_mge.ecsv` OR custom MGE fit | Yes |
| HST F555W image | `data/raw/{galaxy}/hst_*_f555w_*_drz.fits` | For PSF fitting |
| HST F814W image | `data/raw/{galaxy}/hst_*_f814w_*_drz.fits` | For custom MGE |
| SAURON kinematics | `data/raw/{galaxy}/sauron/{galaxy}_r*_idl.fits.gz` | For `_s` files |
| SAURON C2D cube | `data/raw/{galaxy}/sauron/MS_{galaxy}_r*_C2D.fits` | For `_s` files |

### OASIS FITS Format

`read_kinematics_oasis()` in `generate_kin_input.py` expects **HDU 1** with
vector columns (1 row, each cell is a numpy array):

| Field | Shape | Description |
|-------|-------|-------------|
| `XBIN`, `YBIN` | `(n_bins,)` | Bin center coordinates (arcsec) |
| `VEL`, `DVEL` | `(n_bins,)` | Velocity & error |
| `SIG`, `DSIG` | `(n_bins,)` | Velocity dispersion & error |
| `H3`–`H6`, `DH3`–`DH6` | `(n_bins,)` | Gauss-Hermite moments & errors |
| `XPIX`, `YPIX` | `(n_pixels,)` | Pixel coordinates (arcsec) |
| `BINNUM` | `(n_pixels,)` | Pixel-to-bin assignment (0-indexed) |
| `SURF` | `(n_pixels,)` | Surface brightness (can be all-zero) |

The pipeline auto-subtracts systemic velocity via `pafit.fit_kinematic_pa`,
so the FITS **must** contain observed-frame velocities. Pre-subtracted
velocities will still work (the fitted residual will be ~0).

Other formats (CALIFA, ATLAS3D, NIFS) have their own reader functions in
`generate_kin_input.py`. See `kin_input` parameter of `create_kin_input()`.

## Pipeline Overview

```
Step 1: Data Preprocessing
───────────────────────────
raw OASIS FITS.gz                    raw SAURON C2D + PXF_bin FITS
    │                                     │
    └── prep_oasis_batch.py               └── prep_sauron_batch.py
           (scripts/)                           (scripts/)
           aperture_o / bins_o / kin_o         aperture_s / bins_s / kin_s
           mge.ecsv / description.yaml         + sigma ×1.05
    │                                     │
    │                                     │
    └── OASIS + SAURON ── combine ──► gauss_hermite_kins.ecsv

Note: create_kin_input() does NOT write PSF metadata. PSF is injected
separately via add_psf_to_datafile() (see prep_sauron_batch.py Step 3).
The combine step reads PSF from each input ECSV to build the merged metadata.

Step 2: HST + PSF
─────────────────
Download HST F555W ──► fit_psf.py (OASIS)  → psf_fit/psf_result.yaml
                        fit_psf.py (SAURON) → psf_fit_sauron/psf_result.yaml

Step 3: JAM Config + Run
────────────────────────
Create JAM config ──► fill PSF from Step 2 ──► run_jam.py
                                                      │
                                                      └──► results/JAM/{galaxy}/{model}/
```

## Step 1: Data Preprocessing

Decompress OASIS kinematics, generate aperture/bins/kinematics ECSV files, extract
MGE from the SAURON catalog, and write a `description.yaml` stub.

```python
# Single galaxy
import sys; sys.path.insert(0, 'scripts')
from prep_oasis_batch import prep_galaxy
prep_galaxy("NGC4552")

# The script lives in: scripts/prep_oasis_batch.py
```

After this step, `data/processed/{galaxy}/` contains:

| File | Content |
|------|---------|
| `aperture_o.dat` | Bounding box, rotation angle, grid dimensions |
| `bins_o.dat` | Flattened 2D Voronoi bin map |
| `gauss_hermite_kins_o.ecsv` | GH kinematics (v, σ, h3–h6) per bin |
| `mge.ecsv` | MGE photometric decomposition |
| `description.yaml` | Galaxy metadata stub (see below) |

### SAURON preprocessing (`_s` suffix)

The ATLAS3D SAURON data provides wider-field kinematics that can be combined
with OASIS inner kinematics. Requires C2D cube and PXF_bin FITS in
`data/raw/{galaxy}/sauron/`.

```bash
cd scripts
python prep_sauron_batch.py
```

This runs all 3 galaxies (NGC4552/5846/5813) by default. Per galaxy:
- Outputs: `aperture_s.dat`, `bins_s.dat`, `gauss_hermite_kins_s.ecsv`
- SAURON sigma is offset by ×1.05 (known calibration factor)
- PA is automatically checked against OASIS; discrepancy >5° triggers a WARNING
- PSF metadata is auto-read from `results/JAM/{galaxy}/psf_fit_sauron/`
- Combines OASIS (`_o`) and SAURON (`_s`) into `gauss_hermite_kins.ecsv`

### PA Convention

The same angle appears in three places with different representations:

| Location | Convention | Example (NGC5813) |
|----------|-----------|-------------------|
| `pafit` output / `description.yaml` `pa` | Astronomical PA (N→E, degrees) | 147.0° |
| `aperture_o.dat` line 4 | `90° − PA`, CCW from major axis to data X-axis | −57.0° |
| SAURON `angle_deg` input | Same as aperture convention (`90° − pa_o`) | −57.0° |

All are equivalent; `create_kin_input` handles the conversion internally.

### description.yaml format

Auto-generated stub with embedded McDermid+2006 Table 1 data. Fields marked
`TODO` require manual review. Append-only sections (`processing_log`, `jam_models`)
are updated by subsequent pipeline steps.

```yaml
galaxy: NGC4552
type: TODO
distMpc: 15.9
fwhm: 0.67
mge:
  source: sauron_mge.ecsv
  components: 11
  q_min: 0.80
kinematics:
  bins: 709
  pa: 122.5
data_quality: |
  TODO: assess overall quality, note any issues.
processing_log: []
jam_models: []
```

To regenerate just the stub without re-running the full preprocessing:

```bash
python .agents/skills/galaxy-data-prep/scripts/write_description.py NGC4552
```

This script is pipeline-agnostic — it only reads already-generated files
(mge.ecsv, gauss_hermite_kins_o.ecsv, aperture_o.dat) and writes a YAML stub.

### Updating description.yaml

Use `update_description.py` for all subsequent updates. It writes to the YAML
file directly — no manual YAML editing, no temporary scripts. The same module is
also importable from other Python scripts.

```bash
# ---- processing_log entries ----

# Append a log entry
python .agents/skills/galaxy-data-prep/scripts/update_description.py NGC4552 \
    --log "PSF fit" "sigma1=0.291, sigma2=1.320, chi2=0.021"

# With a custom date
python .agents/skills/galaxy-data-prep/scripts/update_description.py NGC4552 \
    --log "manual fix" "replaced HST image" --date 2026-05-01

# ---- jam_models entries ----

# Append a new model
python .agents/skills/galaxy-data-prep/scripts/update_description.py NGC4552 \
    --jam '{"name":"psf_free","chi2":2114,"q":0.695,"lg_mbh":9.09,"lg_ml":0.594,"ratio":1.097,"psf_source":"free"}'

# Replace an existing model (matched by name)
python .agents/skills/galaxy-data-prep/scripts/update_description.py NGC4552 \
    --jam-replace '{"name":"default","chi2":2277,"q":0.515,"lg_mbh":8.80}'

# ---- top-level fields ----

# Set or update any field
python .agents/skills/galaxy-data-prep/scripts/update_description.py NGC4552 \
    --set data_quality "PSF FWHM 0.67. HST image marginal." \
    --set redshift "0.0035" \
    --set type "E0-1"
```

**Python API** (importable from other scripts):

```python
from update_description import append_log, set_field, append_jam_model, replace_jam_model

append_log("NGC4552", step="PSF fit",
           notes="sigma1=0.291, sigma2=1.320, chi2=0.021")

set_field("NGC4552",
          data_quality="HST image marginal, use with caution",
          redshift="0.0035")

append_jam_model("NGC4552",
                 name="psf_free", chi2=2114, q=0.695,
                 lg_mbh=9.09, lg_ml=0.594, ratio=1.097,
                 psf_source="free")
```

**Rule**: when a pipeline step produces an anomaly (large PA discrepancy, failed
PSF fit, suspicious chi2), **ask before updating description.yaml**. Normal
completions can be auto-appended.

## Step 2: HST Image + PSF Fitting

### Download HST F555W Image

Required for PSF fitting. Query the MAST archive:

```python
from astroquery.mast import Observations

# Check availability
obs = Observations.query_criteria(
    target_name="NGC4552", obs_collection="HST",
    filters=["F555W"], dataproduct_type="image"
)

# The combined HLA product filename follows this pattern:
# hst_{zero_padded_PID}_{obs_id}_wfpc2_f555w_pc_drz.fits
```

**Direct download** — the MAST API returns a `mast:` URI, but the actual HTTPS
download URL is:

```
https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLA/url/cgi-bin/getdata.cgi?dataset={filename}
```

Save to `data/raw/{galaxy}/` alongside the OASIS files.

**No HST data?** Some galaxies (e.g., NGC4262) have no HST observations. In that
case, skip PSF fitting and use a rough estimate from the literature (FWHM ÷ 2.355)
or set PSF as a free parameter in JAM.

**MAST returns 503?** The MAST Portal API occasionally returns 503. Wait a few
minutes and retry, or try the direct download URL which bypasses the portal service.

### Fit Multi-Gaussian PSF

Matches HST F555W image to IFU flux map by convolving with a 2-Gaussian PSF.
Uses CMA-ES optimization. Supports both OASIS and SAURON flux images.

```bash
cd JAM

# OASIS PSF
python scripts/fit_psf.py NGC4552 [--niter 50] [--ncpu 100]

# SAURON PSF (auto-applies r<15" circular mask to exclude foreground objects)
python scripts/fit_psf.py NGC4552 --kin-input ATLAS3D

# Disable mask (e.g. for OASIS which doesn't need it)
python scripts/fit_psf.py NGC4552 --mask-radius 0
```

Output directories:
- OASIS: `results/JAM/{galaxy}/psf_fit/`
- SAURON: `results/JAM/{galaxy}/psf_fit_sauron/`

Each output dir contains:
- `psf_result.yaml` — best-fit sigma1, sigma2, weight1, xcenter, ycenter, chi2
- `psf_convergence.png` — chi2 vs generation
- `psf_comparison.png` — IFU vs convolved HST side-by-side

After fitting, the results are auto-appended to `description.yaml`.

**SAURON circular mask**: SAURON's wide field (33"×41") often captures foreground
stars/galaxies. The default r<15" mask excludes these. If the SAURON fit still
hits parameter bounds (sigma2→6.0) or chi2 >100, inspect the HST image quality.

## Step 3: JAM Config + Run Model

**Full config spec & parameter reference → `JAM/AGENTS.md`.**

### Draft the config before PSF fitting

Create `JAM/configs/{galaxy}/default.yaml` from an existing template
(e.g. `NGC4552/default.yaml`). Fill in known values first:

```yaml
name: oasis_test          # any descriptive name
align: cyl                # cyl or sph
logistic: false           # use true for sph models
distMpc: 15.9             # distance from McDermid+2006 or literature
pixsize: 0.27             # OASIS pixel size
aperture_expr: _o         # matches file suffix from Step 1
kinematic_settings:
  use_hist: false         # skip histogram lookup (not generated by default)
  h4_corr: false          # disable h4 correction (see JAM/AGENTS.md § h4_corr)
```

### Critical: q upper bound

The intrinsic flattening parameter `q` appears in the inclination formula as
`sqrt(qmin² - q²)`. If the optimizer explores `q > qmin`, the sqrt becomes
imaginary (NaN) and JAM crashes with:

```
AssertionError: Inclination too low q_lum < ...
```

**Always set `q.hi` strictly less than the MGE's minimum observed `q_obs`.**
Check `description.yaml → mge.q_min` after Step 1.

### After PSF fitting → update config

| Source | When to use |
|--------|------------|
| `fit_psf.py` output (Step 2) | Preferred — fits Gaussian PSF to match HST image |
| McDermid+2006 FWHM ÷ 2.355 | Quick estimate, single-Gaussian |
| Free parameter (`fixed: false`) | Good for initial exploration |

Update `sigma1`, `sigma2`, `weight1` with fitted values from `psf_result.yaml`
and set `fixed: true`.

### Output path convention

```yaml
io_settings:
  input_directory: ../data/processed/{galaxy}
  work_directory: ../results/JAM/{galaxy}/{model_name}
```

### Run JAM

```bash
cd JAM
python scripts/run_jam.py configs/{galaxy}/default.yaml
```

Note: `run_jam.py` hardcodes `sampler = 'cmaes'` (20 generations × 50 population)
regardless of the config's `sampler_type` field. Each run takes 5–7 minutes.

After the run, update `description.yaml`:

```yaml
jam_models:
  - name: default
    chi2: 2127
    q: 0.51
    ratio: 1.06
    lg_mbh: 8.94
    lg_ml: 0.60
    psf_source: fitted
```

## Scripts Reference

| Script | Location | Purpose |
|--------|----------|---------|
| `prep_oasis_batch.py` | `scripts/` | Step 1: OASIS pretreatment for a galaxy |
| `prep_sauron_batch.py` | `scripts/` | SAURON preproc: kin→_s files, sigma offset, PSF, combine |
| `write_description.py` | `.agents/skills/galaxy-data-prep/scripts/` | Regenerate description.yaml stub |
| `update_description.py` | `.agents/skills/galaxy-data-prep/scripts/` | Append/update description.yaml fields |
| `generate_kin_input.py` | `Axi_Schwarzschild/data_prep/` | Low-level: aperture/bins/kin from raw FITS |
| `data_combine.py` | `Axi_Schwarzschild/data_prep/` | MGE creation, PSF injection, aperture merging |
| `fit_psf.py` | `JAM/scripts/` | Step 2: multi-Gaussian PSF fitting |
| `run_jam.py` | `JAM/scripts/` | Step 3: run JAM CMA-ES model |
| `hist_losvd.py` | `Axi_Schwarzschild/data_prep/` | Optional: GH → LOSVD histogram conversion |

## New Galaxy with Custom Kinematics

When adding a galaxy that has **no** entry in `data/raw/` (or you want to
use your own pPXF output instead of the standard McDermid+2006 kinematics):

```python
import sys; sys.path.insert(0, "Axi_Schwarzschild/data_prep")
from generate_kin_input import create_kin_input
from data_combine import create_mge_table
from astropy.io import ascii

galaxy = "NGCxxxx"
out_dir = f"data/processed/{galaxy}/"
custom_fits = "ppxf/output/NGCxxxx/ppxf_kins_NGCxxxx_h6.fits"

# 1. Generate _o files directly from custom FITS
create_kin_input(galaxy, custom_fits, out_dir, expr="_o",
                 kin_input="OASIS", fit_PA=True, plot=True,
                 min_gh_err=0.01)

# 2. Extract MGE from the shared SAURON catalog
mge_all = ascii.read("data/processed/sauron_mge.ecsv")
rows = mge_all[mge_all["galaxy"] == galaxy]
create_mge_table(out_dir, rows["I"], rows["sigma"], rows["q"])

# 3. Write description.yaml stub (manual or via write_description.py)
# python .agents/skills/galaxy-data-prep/scripts/write_description.py {galaxy}
```

Differs from the standard `prep_oasis_batch.py` flow in that it:
- Skips the FITS.gz decompress step (kinematics FITS is already generated)
- Uses any FITS path instead of the hardcoded `data/raw/{galaxy}/oasis/` layout
- Same MGE extraction and `description.yaml` creation as the standard path

After this, add SAURON `_s` files (if available) with PSF metadata and
combine — same as the standard flow.

## Replacing OASIS Kinematics (Variant Directory)

To use custom kinematics (e.g. re-fitted pPXF, different template library) while
keeping the same binning, MGE, and SAURON data — without overwriting the
original processed directory:

```bash
GALAXY="NGCxxxx"
TAG="h6"   # or any suffix you prefer
ORIG="data/processed/$GALAXY"
NEW="data/processed/$GALAXY-$TAG"

# 1. Create variant directory, copy shared files
mkdir -p "$NEW"
cp "$ORIG"/mge.ecsv "$NEW/"
cp "$ORIG"/bins_s.dat "$ORIG"/aperture_s.dat "$ORIG"/gauss_hermite_kins_s.ecsv "$NEW/"
```

```python
# 2. Generate new _o files from custom FITS
import sys; sys.path.insert(0, 'Axi_Schwarzschild/data_prep')
from generate_kin_input import create_kin_input

GALAXY = "NGCxxxx"
TAG = "h6"
create_kin_input(
    galaxy=GALAXY,
    file="path/to/your_custom_ppxf.fits",
    dyn_model_dir=f"data/processed/{GALAXY}-{TAG}/",
    expr="_o", kin_input="OASIS",
    fit_PA=True, plot=True, min_gh_err=0.01,
)
```

```python
# 3. Add PSF metadata & combine with SAURON
from data_combine import add_psf_to_datafile, combine_kin_file

OUT = f"data/processed/{GALAXY}-{TAG}"
# Read PSF from results/JAM/{GALAXY}/psf_fit/psf_result.yaml
add_psf_to_datafile(
    sigma=[s1, s2], weight=[w1, 1 - w1],
    datafile=f"{OUT}/gauss_hermite_kins_o.ecsv",
)
combine_kin_file(f"{OUT}/", ["_o", "_s"])
```

The original `data/processed/{GALAXY}/` remains untouched. The new directory
can be used directly in JAM configs via `io_settings.input_directory`.

## Common Issues

**`AssertionError: Inclination too low`** — q upper bound ≥ MGE qmin. Reduce
`q.hi` in the config to be strictly less than `mge.q_min` from `description.yaml`.

**`KeyError: 'lg_ml_inf'`** — The JAM PARAMS list requires this field even when
`varying_ml` is not used. Add it to the config with `fixed: true`.

**`AttributeError: 'PosixPath' object has no attribute 'find'`** — `use_hist` is
`true` but no histogram file exists. Set `kinematic_settings.use_hist: false`.

**HST image data quality** — Some galaxies have problematic HST images
(e.g. NGC4552: F555W 5a variant has 42% negative pixels). Both OASIS and
SAURON PSF fits may be unreliable. JAM dynamics can shift significantly
with PSF freed. After downloading an HST image, inspect it for negative
pixel fractions and any visible artifacts before PSF fitting.

**SAURON PA mismatch** — OASIS and SAURON velocity field PAs may differ at
large radii. `create_kin_input` auto-checks and prints the difference. A
discrepancy >30° (e.g., NGC5846, 59°) indicates strong kinematic twist —
the galaxy is NOT axisymmetric and should be excluded from axisymmetric
JAM/AxiSchw models.

**SAURON PSF fit fails without mask** — The wide SAURON field (33"×41") often
includes foreground stars/galaxies. `fit_psf.py --kin-input ATLAS3D` defaults
to r<15" circular mask. If chi2 remains >100, inspect the HST image quality.

**MAST 503 errors** — The Portal API's `/columnsconfig` endpoint is unreliable.
Use the direct download URL (`mast.stsci.edu/api/v0.1/Download/file?...`) which
bypasses the portal service.

**No HST data for a galaxy** — Query with `astroquery.mast.Observations` first.
If no observations exist (e.g., NGC4262, NGC4564), skip Steps 3–4 and set PSF
as a free parameter in the JAM config, or estimate from FWHM ÷ 2.355.
