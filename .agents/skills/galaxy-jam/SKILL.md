---
name: galaxy-jam
description: >
  Create and run JAM (Jeans Anisotropic MGE) models for axisymmetric galaxies
  using CMA-ES optimization. Use when fitting JAM dynamical models to estimate
  black hole mass, M/L, inclination, and anisotropy from IFU kinematics. Covers
  config creation, execution, and result logging. Requires preprocessed galaxy
  data from the galaxy-data-prep skill.
---

# Galaxy JAM Modeling

Creates and runs JAM (Jeans Anisotropic MGE) models for axisymmetric galaxies
using CMA-ES optimization.

## Prerequisites

- Preprocessed galaxy data in `data/processed/{galaxy}/` (see `galaxy-data-prep` skill)
- PSF fitting results in `results/JAM/{galaxy}/psf_fit/`
- Conda env: `schw`

## Config & Run

### Create config

Write `JAM/configs/{galaxy}/default.yaml` from the NGC5813-h6 template:

```yaml
name: default
align: cyl                    # cyl or sph
logistic: false
nstep: 20000
distMpc: 25.0                 # from description.yaml
pixsize: 0.27                 # OASIS pixel size (arcsec)
aperture_expr: _o             # match output from preprocessing
kinematic_settings:
    use_hist: false
    h4_corr: false

parameters:
    q:                         # intrinsic flattening
        value: 0.5
        fixed: false
        grid_setting: {lo: 0.1, hi: 0.25, sigpar: 0.05}
    ratio:                     # velocity anisotropy
        value: 1.0
        fixed: false
        grid_setting: {lo: 0.5, hi: 1.5, sigpar: 0.05}
    lg_mbh:                    # log10(M_BH / M_sun)
        value: 8.6
        fixed: false
        grid_setting: {lo: 6.14, hi: 10.14, sigpar: 0.01}
    lg_ml:                     # log10(M/L)
        value: 0.5
        fixed: false
        grid_setting: {lo: -0.5, hi: 0.7, sigpar: 0.05}
    f_dm:                      # dark matter fraction
        value: 0.0; fixed: true
    lg_rb:                     # NFW break radius
        value: 3.0; fixed: true
    sigma1/sigma2/weight1:     # PSF from Step 2
        value: <from psf_result.yaml>; fixed: true
    ratio_inf/lg_ra/lg_al:     # logistic anisotropy (fixed)
        value: 1.0/0.0/0.0; fixed: true
    lg_ml_inf:                 # varying M/L (fixed)
        value: 0.5; fixed: true

parameter_space_settings:
    sampler_type: Adamet       # ignored; run_jam.py hardcodes cmaes
    sampler_settings:
        seed: 42
        n_processes: 16
io_settings:
    input_directory: ../data/processed/{galaxy}
    work_directory: ../results/JAM/{galaxy}/{model_name}
```

### Sph logistic template (sph2 style)

For spherical alignment with logistic radial anisotropy profile (`ratio>1`, `ratio_inf<1`):

```yaml
name: sph2
align: sph
logistic: true
nstep: 20000
distMpc: 25.0
pixsize: 0.27
aperture_expr: _o
kinematic_settings:
    use_hist: false
    h4_corr: false

parameters:
    q:
        value: 0.2; fixed: false
        grid_setting: {lo: 0.1, hi: 0.25, sigpar: 0.05}
    ratio:                     # σ_θ/σ_r > 1 → tangentially biased at center
        value: 1.1; fixed: false
        grid_setting: {lo: 1.0, hi: 2.0, sigpar: 0.1}
    lg_mbh:
        value: 8.6; fixed: false
        grid_setting: {lo: 8.0, hi: 10.0, sigpar: 0.01}
    lg_ml:
        value: 0.5; fixed: false
        grid_setting: {lo: -0.5, hi: 0.7, sigpar: 0.05}
    f_dm:
        value: 0.0; fixed: true
    lg_rb:
        value: 3.0; fixed: true
    sigma1/sigma2/weight1:
        value: <from psf_result.yaml>; fixed: true
    ratio_inf:                 # σ_θ/σ_r < 1 → less tangential / radial at ∞
        value: 0.9; fixed: false
        grid_setting: {lo: 0.5, hi: 1.0, sigpar: 0.1}
    lg_ra:
        value: 0.0; fixed: true
    lg_al:
        value: 0.3; fixed: true
    lg_ml_inf:
        value: 0.5; fixed: true

parameter_space_settings:
    sampler_type: Adamet
    sampler_settings: {seed: 42, n_processes: 16}

io_settings:
    input_directory: ../data/processed/{galaxy}
    work_directory: ../results/JAM/{galaxy}/sph2
```

5 free parameters: q, ratio, lg_mbh, lg_ml, ratio_inf.

### Critical: q upper bound

The intrinsic flattening `q` appears in the inclination formula as
`sqrt(qmin² - q²)`. If `q > qmin`, the sqrt becomes imaginary and JAM
crashes:

```
AssertionError: Inclination too low q_lum < ...
```

**Always set `q.hi` strictly less than `mge.q_min`** from `description.yaml`.
E.g. if MGE `q_min = 0.263`, set `q.hi = 0.25`.

### PSF values

Update `sigma1`, `sigma2`, `weight1` from `fit_psf.py` output
(`psf_result.yaml`). Leave them `fixed: true`.

### Run JAM

```bash
# Always use setsid to fully detach from shell
cd JAM && setsid /home/haitong/miniconda3/envs/schw/bin/python -u \
    scripts/run_jam.py configs/{galaxy}/default.yaml \
    > /tmp/{galaxy}_jam.log 2>&1 &
```

**Note:** `nohup` alone is insufficient — the process may be killed when the
parent shell exits. Use `setsid` instead.

Run time depends on alignment:
- **cyl** (analytic LOS integration): ~5–7 min
- **sph** (numerical LOS integration): ~25–30 min

### After the run

#### Compute inclination

```python
import math
q_min = 0.263   # from description.yaml mge.q_min
q = 0.219       # from best-fit
inc = math.degrees(math.asin(math.sqrt((1 - q_min**2) / (1 - q**2))))
```

#### Update description.yaml

```yaml
jam_models:
  - name: default
    chi2: 490.2
    q: 0.100
    lg_mbh: 8.89
    lg_ml: 0.588
    ratio: 0.956
    inc: 75.8
    psf_source: fitted
    kinematic_source: pPXF h6 → 0±0.1 prior
  - name: sph2
    chi2: 392.2
    q: 0.219
    lg_mbh: 8.96
    lg_ml: 0.550
    ratio: 1.404
    ratio_inf: 0.720
    inc: 81.4
    psf_source: fitted (proposal 6099)
    kinematic_source: pPXF h6 → 0±0.1 prior
```

Include `ratio_inf` for logistic models and `inc` for all models. Use `update_description.py` or edit the YAML directly.

修改后提交到根仓库：

```bash
git add data/processed/{galaxy}/description.yaml
git commit -m "{galaxy}: add {model} results"
git push
```

### Cleaning __pycache__

If you edit `jam_fit/jam.py` or other JAM source files, stale `.pyc` can
cause mysterious errors (e.g. wrong function signatures). Clear caches:

```bash
rm -rf JAM/jam_fit/__pycache__ JAM/__pycache__
```

## Scripts Reference

| Script | Location | Purpose |
|--------|----------|---------|
| `run_jam.py` | `JAM/scripts/` | JAM CMA-ES model runner |

Full config spec → `JAM/AGENTS.md`.
