---
name: dyn-workflow
description: >
  Full Dynamite triaxial Schwarzschild workflow: create configs, submit model
  runs to Slurm, monitor progress, recover from failures, and post-process
  results. Use when starting a new triaxial model run or analyzing completed
  Dynamite model grids.
---

# Dynamite (Trischwarzpy) Workflow

## Quick start

```bash
# Submit a new run
export PATH=/soft/slurm/bin:$PATH
cd Trischwarzpy
python scripts/run_bo_dynamite.py ../dyn_config/<galaxy>.yaml
```

---

## Step 0: Create config & submit

### 0a. Create config

Configs are Dynamite-native YAML with `system_components`. Use an existing config as template:

```bash
cp ../dyn_config/5813-hist-bayes.yaml ../dyn_config/<new>.yaml
```

Required edits: `name`, `input_directory` (templates/), `output_directory` (dyn_models/), parameters.

Key config sections:

- `system_attributes`: `distMPc`, `name`, `position_angle`
- `system_components.bh`: Plummer BH, `m` (log M_BH), `a` (fixed 0.001)
- `system_components.dh`: NFW dark halo, `include: false` to disable
- `system_components.stars`: TriaxialVisibleComponent, `mge_pot`, `mge_lum`, `q/p/u` shape params
- `system_components.stars.kinematics`: {oasis, atlas3d} with type `GaussHermite` or `BayesLOSVD`
- `orblib_settings`: `nE`, `nI2`, `nI3`, `dithering`, `logrmin`, `logrmax`
- `weight_solver_settings`: `NNLS` with `nnls_solver: scipy`, `regularisation`, `number_GH`
- `parameter_space_settings.generator_type`: `BayesOpt` or `SpecificModels`
- `parameter_space_settings.generator_settings`: `n_processes`, `strategy` (`cl_mean`|`gp`)
- `parameter_space_settings.stopping_criteria`: `n_max_mods`, `n_max_iter`
- `multiprocessing_settings`: `ncpus` (pool workers), `ncpus_weights` (concurrent NNLS)
- `io_settings`: `input_directory` → `../templates/<galaxy>`, `output_directory` → `../dyn_models/<model>`

### 0b. Submit to Slurm

```bash
export PATH=/soft/slurm/bin:$PATH
cd /share/home/maoshudeLab/wanght245001/galactic_dynamics/Trischwarzpy
sbatch submit/<script>.sh
```

SLURM template:
```bash
#!/bin/bash
#SBATCH -p test,test-intel
#SBATCH -N 1 -n 1 -c 32
#SBATCH --mem=300g
cd $SLURM_SUBMIT_DIR
source ~/miniconda3/etc/profile.d/conda.sh
conda activate schw
python scripts/run_bo_dynamite.py ../dyn_config/<name>.yaml
```

### 0c. Resume from existing output

If the job crashed mid-run, kill it and resume without losing completed models:

```bash
scancel <jobid>
# Fix any issues (config, resources, etc.)
sbatch submit/<script>.sh    # script uses -r flag
```

The `-r`/`--resume` flag sets `reset_existing_output=False`, preserving `all_models.ecsv` and continuing from the last completed iteration.

### 0d. Monitor progress

```bash
# Job status
squeue -u wanght245001

# Per-iteration progress (best chi2, top models)
SCHW=/share/home/maoshudeLab/wanght245001/miniconda3/envs/schw/bin/python3
$SCHW Trischwarzpy/scripts/check_progress.py dyn_models/<name>

# All models at once
$SCHW Trischwarzpy/scripts/check_progress.py --all

# Continuous watch (every 30s)
$SCHW Trischwarzpy/scripts/check_progress.py --all --watch 30

# Watch log tail
tail -f Trischwarzpy/log/<name>.err
```

### Common issues

| Symptom | Cause | Fix |
|---------|-------|-----|
| `YAML ParserError` at `\mathrm` | LaTeX braces in flow-style `{...}` break YAML parser | Use block-style YAML (multi-line keys), never inline `{key: val}` format |
| Job runs 3 seconds, exits | Config has `reset_existing_output=True` and output dir exists | Use `-r` flag to resume, or delete output dir |
| `ModuleNotFoundError: No module named 'pathos'` | Using base `python3` instead of `schw` env python | Always use `/share/.../miniconda3/envs/schw/bin/python3` or `conda activate schw` first |
| Pool hangs on last few models | Worker OOM-killed (NNLS memory blowup) leaks NNLS semaphore | Kill job, keep `ncpus=64`, reduce `ncpus_weights=4`, resume with `-r` |
| OOM kills in SLURM log: `Detected N oom_kill events` | Concurrent NNLS solves exceed memory | Each NNLS uses **3-5× histogram file size** (~40GB per 9GB histogram). `ncpus_weights=8` → 320GB peak, exceeds 300GB. **Set `ncpus_weights=4`** (160GB safe). `ncpus=64` is fine (orbit integration is memory-light). |
| NaN values in all_models.chi2 | `orblib_done=False` or model crashed | Run recovery script (Step 1) |
| `No previous models found` | `reset_existing_output=True` with empty output dir | Normal for first run |
| GP contour delta χ² too small (0-5 range) | chi2 passed as `np.log10(chi2)` instead of real chi2 | Use real `kinchi2` values. Log scaling hides the parameter landscape and makes 1σ/2σ/3σ contour lines meaningless. |
| `python3 -c` commands time out in `ssh_helper.py` | Base python, not schw env, missing `astropy`/`pathos` | Use `conda activate schw && python3 -c "..."` or full path to schw python |

---

## Step 1: Recover models after crash

When the pool hangs with completed models on disk but `all_models.ecsv` out of date:

```python
import os, numpy as np
from astropy.table import Table

home = '../dyn_models/<name>'
t = Table.read(f'{home}/all_models.ecsv', format='ascii.ecsv')
n = 0
for row in t:
    if row['all_done']:
        continue
    wf = f"{home}/models/{row['directory']}/orbit_weights.ecsv"
    if os.path.exists(wf):
        w = Table.read(wf, format='ascii.ecsv')
        row['chi2'] = w.meta.get('chi2_tot', np.nan)
        row['kinchi2'] = w.meta.get('chi2_kin', np.nan)
        row['orblib_done'] = True
        row['weights_done'] = True
        row['all_done'] = True
        n += 1
t.write(f'{home}/all_models.ecsv', format='ascii.ecsv', overwrite=True)
print(f'Recovered {n} models')
```

Then resume with `-r`.

---

## Step 2: Post-processing

### 2a. Load results

```python
import numpy as np
from astropy.table import Table
from Trischwarzpy.mod_dyn import Configuration
import dynamite as dyn

# Load config
c = Configuration('../dyn_config/<name>.yaml', reset_logging=False, reset_existing_output=False)
plotter = dyn.plotter.Plotter(config=c)

# Table of all models
t = plotter.all_models.table
done = t[t['all_done'] == True]

# Best models
idx = plotter.all_models.get_best_n_models_idx(n=10, which_chi2='kinchi2')
for i in idx:
    row = t[i]
    print(f'bh={row["m-bh"]:.2e}  ml={row["ml"]:.2f}  kinchi2={row["kinchi2"]:.0f}')
```

### 2b. Chi2 analysis

**Use `analyze_results.py` (recommended):**

```bash
# Must use schw env python (not base python3) — pathos is installed there
SCHW=/share/home/maoshudeLab/wanght245001/miniconda3/envs/schw/bin/python3
$SCHW Trischwarzpy/scripts/analyze_results.py dyn_models/<name> -o dyn_models/<name> --gp
```

Generates: `chi2_landscape.png`, `convergence.png`, `histograms.png`, `gp_contour.png`.

**GP contour manually:**

```python
from Trischwarzpy.mod_dyn.postproc import plot_gpcontour

pars = np.log10(np.vstack([done['m-bh'], done['ml']]).T)
chi2 = done['kinchi2']        # use REAL chi2, NOT log10!
parnames = [r'$\log M_{\rm BH}$', r'$\log M/L$']
fig, axes = plot_gpcontour(pars, chi2, parnames=parnames, levels=20)
```

> **Important**: Pass real `kinchi2` values, not `np.log10(kinchi2)`. Log scaling makes the delta χ² scale meaningless and the 1σ/2σ/3σ black dashed contour lines (at Δχ²=2.3/6.17/11.8) won't work.

**GP contour with Dynesty integration (uncertainty estimation):**

```python
from Trischwarzpy.mod_dyn.postproc import gpfit, dynesty_sample

pars = np.log10(np.vstack([done['m-bh'], done['ml']]).T)
chi2 = done['kinchi2']
gp = gpfit(normalize_parlist(pars), chi2, noise=30)
par_range = np.percentile(pars, [0, 100], axis=0).T

# use_log=False means chi2 is real (not log-transformed)
results, rlist = dynesty_sample(gp, par_range, gpsamples=10, nlive=1000, use_log=False)

from dynesty import plotting as dyplot
cfig, caxes = dyplot.cornerplot(
    results, labels=parnames,
    quantiles=(0.16, 0.5, 0.84), show_titles=True,
    title_fmt='.3f', quantiles_2d=1.0 - np.exp(-0.5 * np.arange(3.1)**2))
```

### 2c. Re-run weight solver (redo_weight.py)

Useful for trying different NNLS solvers or re-solving failed models:

```bash
# Edit redo_weight.py to specify model range, then:
python scripts/redo_weight.py ../dyn_config/<name>.yaml
```

Or programmatically:

```python
from Trischwarzpy.mod_dyn.helper import resolve_weight, save_kinmaps_data, update_chi2_table
import dynamite as dyn

c = Configuration('../dyn_config/<name>.yaml', reset_logging=False, reset_existing_output=False)
plotter = dyn.plotter.Plotter(config=c)
m = plotter.all_models.get_model_from_row(idx)  # dyn.model.Model

resolve_weight(m, solver='scipy')                    # re-do NNLS
save_kinmaps_data(c, m, n_aper=2)                   # re-generate kinematic maps
update_chi2_table('../dyn_models/<name>', expr='-new')  # update all_models.ecsv
```

### 2d. Kinematic maps

```python
# Built-in (via plotter, generates maps for all models)
plotter.plot_kinematic_maps(kin_set='all', cbar_lims='default')

# For a specific model
m = plotter.all_models.get_model_from_row(idx)
# .kmod holds the fitted kinematics arrays
```

### 2e. Anisotropy

```python
from Trischwarzpy.mod_dyn.helper import anisotropy_single

m = plotter.all_models.get_model_from_row(idx)
sigmas, moments, rr = anisotropy_single(m)  # returns spherical velocity dispersions
```

### 2f. Compare axisymmetric vs triaxial results

```python
# Load AxiSchw model
from schwarzpy.model import Model, read_ap, plot_vel_ellipsoid
m_axi = Model(axi_work_dir, aperture_exprs=['_o', '_s'])

# Load Dynamite/triaxial model
m_dyn = plotter.all_models.get_model_from_row(idx)

# Compare kinematics
from Trischwarzpy.mod_dyn.helper import save_mer_f
def plot_kincompare(kins1, kins2, labels=['Axisym', 'Triaxi'], nrow=3, ncol=2):
    fig, axes = plt.subplots(nrow, ncol, figsize=(2*ncol*2,2*nrow*2))
    for igh in range(nrow*ncol):
        plt.sca(axes[igh//ncol, igh%ncol])
        plt.scatter(kins1[:,igh], kins2[:,igh], s=5, alpha=0.5)
        plt.plot([kins1.min(), kins1.max()], [kins1.min(), kins1.max()], 'r--')
        label = ['V','sig','h3','h4','h5','h6'][igh]
        plt.title(label); plt.xlabel(labels[0]); plt.ylabel(labels[1]); plt.grid()
```

---

## Quick reference

| Task | Command / Code |
|------|---------------|
| Submit new run | `cd Trischwarzpy && python scripts/run_bo_dynamite.py ../dyn_config/<name>.yaml` |
| Submit with SLURM | `sbatch submit/<script>.sh` |
| Resume from crash | `python scripts/run_bo_dynamite.py -r ../dyn_config/<name>.yaml` |
| Check progress | `$SCHW Trischwarzpy/scripts/check_progress.py --all` |
| Watch continuously | `$SCHW Trischwarzpy/scripts/check_progress.py --all --watch 30` |
| Monitor OOM | `grep oom_kill Trischwarzpy/log/<name>.err` |
| Recover models (post-crash) | Run recovery script (Step 1), then resume with `-r` |
| Analyze results | `$SCHW Trischwarzpy/scripts/analyze_results.py dyn_models/<name> -o dyn_models/<name> --gp` |
| Kinematic maps | `$SCHW Trischwarzpy/scripts/analyze_results.py dyn_models/<name> --kinmap --axi galaxy_models/<axi_name>` |
| Load results | `Configuration(config, reset_existing_output=False)` + `dyn.plotter.Plotter(config=c)` |
| Best model | `plotter.all_models.get_best_n_models_idx(n=10)` |
| GP contour | `plot_gpcontour(pars, chi2, parnames=...)` — use real chi2, not log10! |
| Dynesty sampling | `dynesty_sample(gp, par_range, gpsamples=10, nlive=1000, use_log=False)` |
| Resolve weights | `resolve_weight(model, solver='scipy')` |
| Update chi2 table | `update_chi2_table(home_dir, expr='-new')` |
| Kinematics maps | `plotter.plot_kinematic_maps(kin_set='all')` |
| Anisotropy | `anisotropy_single(model)` → `(sigmas, moments, rr)` |
| OOM-safe config | `ncpus: 64, ncpus_weights: 4` (orbit is memory-light, NNLS needs 3-5× histogram RAM) |
| Compare Axi vs Triax | `plot_kincompare(m_axi.kmod, m_dyn.kmod, ['Axisym','Triaxi'])` |
