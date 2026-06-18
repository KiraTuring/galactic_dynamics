---
name: schw-post-process
description: >
  Post-process AxiSchw Schwarzschild grid results. Use when inspecting completed
  model grids, analyzing chi2 convergence, plotting kinematics comparison
  (V/sigma/h3/h4 Model vs Observed), diagnosing mass & anisotropy profiles,
  running isotropy sensitivity tests with rerun_nnls.py, or evaluating a
  BayesOpt/CmaEs run. Do NOT use for data preparation (use galaxy-data-prep
  instead) or for submitting new model runs.
---

# Schwarzschild Post-Processing

Analyze completed AxiSchw grid optimization outputs.

## When to use

- You have a **finished grid run** (BayesOpt / CmaEs / NestedSamp) and need to
  inspect results.
- User asks about **chi2 convergence**, **best-fit parameters**, or **kinematics
  comparison**.
- User wants to **plot kinematics maps**, **anisotropy profiles**, or **chi2
  contours**.
- User wants to **rerun the weight solver** with different isotropy or mass
  constraints.

## When NOT to use

- Preparing input data → use `galaxy-data-prep`
- Submitting new model runs → use `scripts/run.py` or Slurm directly via AGENTS.md
- Checking pipeline architecture → read AGENTS.md

## Prerequisites

All commands MUST run from the Axi_Schwarzschild repo root — relative paths
in config files are resolved against CWD, not the config location.

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate schw
cd /share/home/maoshudeLab/wanght245001/galactic_dynamics/Axi_Schwarzschild
```

### Path rules

- `Iterator()`: pass the **home directory** (relative or absolute)
- `Model()`: pass the **full work_dir path** (the `bh*_ml*_...` subdirectory)
- `os.chdir(repo_root)` before any `from schwarzpy import ...`

---

## Step 1: Load results

```python
import os, numpy as np
os.chdir('/path/to/Axi_Schwarzschild')  # MUST be repo root!
from schwarzpy.grid import Iterator

# Grid results: pass home directory
iter0 = Iterator('../galaxy_models/<model_dir>')
model_list = iter0.get_model_list()    # structured array: 'dir', 'par', 'chi2'
chi2list = np.array(model_list['chi2'], dtype=float)   # NOTE: chi2_kin, NOT chi2_tot
parlist  = np.array([list(p) for p in model_list['par']], dtype=float)
```

```python
# Single model (multi-aperture: MUST pass aperture_exprs)
from schwarzpy.model import Model
m = Model(work_dir, aperture_exprs=['_o', '_s'])  # OASIS + SAURON (2-aperture)
m = Model(work_dir)                                 # single-aperture
```

> **Pitfalls:**
> - `get_best_model()` creates a Model internally without `aperture_exprs` — fails
>   on multi-aperture. Use `model_list['dir'][0]` instead.
> - `model_list['chi2']` is `chi2_kin` (kinematic only), not `chi2_tot` (kin+mass).
> - Index with `model_list['par'][i]`, not `model_list[i]['par']`.

---

## Step 2: chi2 analysis

```python
# Built-in GP contour plot (returns (fig, axes) tuple)
fig, _ = iter0.plot_grid_chi2(chi2_range=(-10, 500), size=12)
fig.savefig('chi2_landscape.png', dpi=150, bbox_inches='tight')
```

```python
# Alternatively, use postproc for more control
from schwarzpy.postproc import plot_gpcontour
labels = [r'$\log M_{\rm BH}$', r'$\log M/L$']  # adjust per model params
plot_gpcontour(parlist, chi2list, parnames=labels)
```

**Convergence check**: look at `chi2plot-*.png` files in the model home directory.
If chi2 still dropping in the last iteration, increase `max_iterations` in config.

---

## Step 3: Kinematics comparison

Generate the kinematics output file first (one-time cost per model set, does NOT
rebuild orbits):

```bash
python scripts/rerun_nnls.py -n 8 galaxy_models/<model_dir>/bh*/
```

Then plot:

```python
m = Model(work_dir, aperture_exprs=['_o', '_s'])

# Standard panels: V / sigma / h3 / h4 — Model vs Observed vs Residual
m.plot_kinematics()

# 1D major-axis slit extraction (angle in degrees)
m.plot_kin1d(angle=90, slit=2)
```

**Output check**: the residuals panel should show near-zero, symmetric scatter.
Large systematic residuals indicate poor fit or wrong systemic parameters.

---

## Step 4: Mass & anisotropy

```python
# Radial velocity anisotropy beta(r)
m.plot_anisotropy()

# Meridional-plane orbital moments for velocity ellipsoids
from schwarzpy.model import read_nn_mer, read_ap
mer = read_nn_mer(work_dir)
```

---

## Step 5: Isotropy sensitivity

Rerun the weight solver with isotropy constraints to test solution robustness
against this prior.

```bash
# Single model with isotropy
python scripts/rerun_nnls.py -i 0.02 <model_dir>

# Batch (8 parallel), different iso strengths
python scripts/rerun_nnls.py -i 0.02 -n 8 galaxy_models/<model_dir>/bh*/
python scripts/rerun_nnls.py -i 0.1  -n 8 galaxy_models/<model_dir>/bh*/

# Combine isotropy with mass-constraint overrides
python scripts/rerun_nnls.py -i 0.02 \
    -c mer_plane_masses_err=0.01 \
    -c proj_plane_masses_err=0.02 \
    -n 8 galaxy_models/<model_dir>/bh*/
```

Each run produces `weights_<suffix>.ecsv` + `datfil_<suffix>/` per model.
Original files are untouched.

Compare chi2 with vs without isotropy. A small increase suggests the solution
is physically plausible; a large jump (>few %) suggests tension with isotropy.

---

## Quick reference

| Task | Command / Code |
|------|---------------|
| Load grid | `Iterator(home_dir).get_model_list()` → `['dir','par','chi2']` |
| Best model dir | `model_list['dir'][0]` (first entry sorted by chi2) |
| Open model | `Model(dir, aperture_exprs=['_o','_s'])` |
| chi2 GP contour | `iter0.plot_chi2_grid()` → returns `(fig, axes)` tuple |
| chi2 raw scatter | `plot_gpcontour(parlist, chi2list, parnames=...)` |
| Kinematics | `Model.plot_kinematics()` |
| 1D slit | `Model.plot_kin1d(angle=90)` |
| Anisotropy | `Model.plot_anisotropy()` |
| Rerun NNLS (single) | `rerun_nnls.py -i 0.02 model_dir/` |
| Rerun NNLS (batch) | `rerun_nnls.py -i 0.02 -n 8 model_dir/*/` |
| Submit to Slurm | `sbatch --wrap="source ... && conda activate schw && python scripts/rerun_nnls.py ..."` |
