---
name: schw-workflow
description: >
  Full AxiSchw Schwarzschild workflow: create configs, submit model runs to Slurm,
  monitor progress, and post-process results. Covers the complete lifecycle from
  config creation through chi2 analysis, kinematics comparison, mass/anisotropy
  diagnostics, and isotropy sensitivity tests. Use when starting a new model run
  or analyzing completed grid results.
---

# AxiSchw Workflow

Complete lifecycle from submission to analysis.

## Quick start

```bash
# Submit a new run
export PATH=/soft/slurm/bin:$PATH
cd Axi_Schwarzschild
python scripts/run.py ../config/<galaxy>.yaml

# Analyze completed results
python scripts/analyze_results.py ../galaxy_models/<model_dir>
```

---

## Step 0: Create config & submit

### 0a. Create config

Copy an existing config for the same galaxy and modify 4 fields:

```bash
cp config/<base>.yaml config/<new>.yaml
# Edit these fields:
#   name:               unique run name
#   template_directory: ../templates/<template>
#   home_directory:     ../galaxy_models/<model_dir>
#   seed:               different integer (to avoid duplicate sampling)
```

### 0b. Submit to Slurm

```bash
export PATH=/soft/slurm/bin:$PATH
cd /share/home/maoshudeLab/wanght245001/galactic_dynamics/Axi_Schwarzschild
python scripts/run.py ../config/<name>.yaml
```

`run.py` submits both an iterator (BayesOpt sampler) and a worker pool (MPI, 64
processes). Workers process models from `queue.txt`; the iterator generates new
parameters each iteration.

### 0c. Monitor progress

```bash
# Job status
squeue -u wanght245001

# Models completed (growing = running)
ls galaxy_models/<home>/ | wc -l

# Queue status (Q=queued, R=running, F=finished)
tail galaxy_models/<home>/queue.txt
```

### Common issues

| Symptom | Cause | Fix |
|---------|-------|-----|
| `QOSMaxCpuPerUserLimit` | CPU quota full from other jobs | Wait, or `scancel` lower-priority jobs |
| `No such file: config/...` | Path relative to Axi_Schwarzschild | Use `../config/<name>.yaml` |
| `MGE file ... not found` | Wrong `template_directory` | Check templates dir exists at that path |
| Workers idle after submit | Pool still queuing | Wait for pool to show `R` in squeue |

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
m = Model(work_dir, aperture_exprs=['_o', '_s'])
```

> **Pitfalls:**
> - `get_best_model()` fails on multi-aperture. Use `model_list['dir'][0]`.
> - `model_list['chi2']` is `chi2_kin`, not `chi2_tot`.

---

## Step 2: chi2 analysis

```python
# Built-in GP contour plot (returns (fig, axes) tuple)
fig, _ = iter0.plot_chi2_grid(size=12)
fig.savefig('chi2_landscape.png', dpi=150, bbox_inches='tight')
```

```python
# Or use postproc for more control
from schwarzpy.postproc import plot_chi2_contour
plot_chi2_contour(parlist, chi2list, parnames=[r'$\log M_{\rm BH}$', r'$\log M/L$'])
```

**Convergence check**: look at `chi2plot-*.png` in the model home directory.
If chi2 still dropping in the last iteration, increase `max_iterations` in config.

---

## Step 3: Kinematics comparison

```bash
# Prerequisite (one-time per model set):
python scripts/rerun_nnls.py -n 8 galaxy_models/<model_dir>/bh*/
```

```python
m = Model(work_dir, aperture_exprs=['_o', '_s'])
m.plot_kinematics()        # V / sigma / h3 / h4: Model vs Observed vs Residual
m.plot_kin1d(angle=90)     # 1D slit extraction
m.plot_2D_anisotropy()     # 2D velocity ellipsoids in meridional plane
```

---

## Step 4: Mass & anisotropy

```python
m.plot_anisotropy()        # beta(r) velocity anisotropy profile
```

---

## Step 5: Isotropy sensitivity

```bash
# Single model
python scripts/rerun_nnls.py -i 0.02 <model_dir>

# Batch (8 parallel), different iso strengths
python scripts/rerun_nnls.py -i 0.02 -n 8 galaxy_models/<model_dir>/bh*/
python scripts/rerun_nnls.py -i 0.1  -n 8 galaxy_models/<model_dir>/bh*/

# With mass-constraint overrides
python scripts/rerun_nnls.py -i 0.02 \
    -c mer_plane_masses_err=0.01 -c proj_plane_masses_err=0.02 \
    -n 8 galaxy_models/<model_dir>/bh*/
```

Output: `weights_<suffix>.ecsv` + `datfil_<suffix>/` per model. Original untouched.

---

## Quick reference

| Task | Command / Code |
|------|---------------|
| Create config | `cp config/base.yaml config/new.yaml` + edit 4 fields |
| Submit run | `run.py ../config/<name>.yaml` |
| Check jobs | `squeue -u wanght245001` |
| Analyze results | `python scripts/analyze_results.py <model_dir>` |
| Load grid | `Iterator(home_dir).get_model_list()` → `['dir','par','chi2']` |
| Best model dir | `model_list['dir'][0]` |
| Open model | `Model(dir, aperture_exprs=['_o','_s'])` |
| chi2 GP contour | `iter0.plot_chi2_grid()` → `(fig, axes)` |
| Kinematics | `Model.plot_kinematics()` |
| 1D slit | `Model.plot_kin1d(angle=90)` |
| 2D ellipsoids | `Model.plot_2D_anisotropy()` |
| Anisotropy | `Model.plot_anisotropy()` |
| Rerun NNLS | `rerun_nnls.py -i 0.02 -n 8 model_dir/*/` |
