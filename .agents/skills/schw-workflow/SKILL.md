---
name: schw-workflow
description: >
  Full AxiSchw Schwarzschild workflow: create configs, submit model runs to Slurm,
  monitor progress, and post-process results. Covers the complete lifecycle from
  config creation through chi2 analysis, kinematics comparison, mass/anisotropy
  diagnostics, and isotropy sensitivity tests. Use when starting a new axisymmetric
  Schwarzschild model run, analyzing completed grid results, or re-running NNLS
  with different isotropy constraints. Make sure to use this skill whenever the
  user mentions AxiSchw, axisymmetric Schwarzschild, or orbit-superposition
  modeling with schwarzpy.
---

# AxiSchw Workflow

Complete lifecycle from submission to analysis.

## Quick start

```bash
# Submit a new run
export PATH=/soft/slurm/bin:$PATH
cd Axi_Schwarzschild
python scripts/run.py ../config/<galaxy>.yaml

# Check status
python scripts/check_run.py ../galaxy_models/<model_dir>

# Analyze completed results
python scripts/analyze_results.py ../galaxy_models/<model_dir> -o results/<name>/

# Download to local for viewing
ls results/<name>/
```

---

## Step 0: Create config & submit

### 0a. Create config

Copy an existing config for the same galaxy and modify:

```bash
cp config/<base>.yaml config/<new>.yaml
```

| Field | How to set |
|-------|-----------|
| `name` | unique run name |
| `template_directory` | `../templates/<template>` |
| `home_directory` | `../galaxy_models/<model_dir>` (check if exists → ask not overwrite) |
| `distMpc` | from description.yaml |
| `BH` value, lo, hi | from JAM best-fit or previous run; `lo` must be below expected minimum |
| `M/L` value, lo, hi | same, ensure range covers JAM best |
| `inc.value`, `inc.fixed` | **check JAM results** for plausible inclination; do not blindly copy 75° |
| `seed` | different integer per run (avoid duplicate sampling) |
| `clear_existing_models` | `true` for new run, `false` for resume |

> **⚠️ MUST confirm with user before submitting** — present the config parameters
> (name, distMpc, BH/ML ranges, inc, output dir) and wait for explicit approval.
> AxiSchw runs take hours to days — accidental overwrite or wrong parameters
> waste significant time and compute resources.

### 0b. Submit to Slurm

```bash
export PATH=/soft/slurm/bin:$PATH
cd /share/home/maoshudeLab/wanght245001/galactic_dynamics/Axi_Schwarzschild
python scripts/run.py ../config/<name>.yaml
```

`run.py` submits both an iterator (BayesOpt sampler) and a 64-process MPI pool.
Workers process models from `queue.txt`; the iterator generates new parameters
each iteration.

### 0c. Monitor progress

```bash
# Quick status
python scripts/check_run.py ../galaxy_models/<model_dir>

# Job status
squeue -u wanght245001

# Iteration progress
ls galaxy_models/<home>/chi2plot-*.png
```

### MPI pool fails: fall back to independent workers

If the MPI pool doesn't start (PD with `QOSMaxCpuPerUserLimit`), **ask the user**
before submitting independent workers — they consume job slots and may exceed
the QOS job limit. If approved:

```bash
Q="../galaxy_models/<home>/queue.txt"
for i in $(seq 0 15); do
    sbatch -p test --cpus-per-task 1 --mem 8G -J "w-$i" \
        --wrap "source ~/miniconda3/etc/profile.d/conda.sh && conda activate schw && python scripts/schw_proc.py -p $i -q $Q"
done
```

### Common issues

| Symptom | Cause | Fix |
|---------|-------|-----|
| `QOSMaxCpuPerUserLimit` (PD) | 64-CPU MPI pool exceeds user quota | Cancel lower-priority jobs, or use independent workers above |
| `mpirun: command not found` | `module` not sourced in non-login shell | Code now auto-sources; verify with `srun -p test -n 2 --time=2 bash -c "source /etc/profile.d/module-profile.sh && module load openmpi/4.1.8 && which mpirun"` |
| `module: command not found` | `/etc/profile.d/module-profile.sh` not sourced | `source /etc/profile.d/module-profile.sh` before `module load` |
| `Connection closed by UNKNOWN port` | EasyConnect VPN session expired | Re-login at http://localhost:8080 (VNC password: `opencode`) |
| `No such file: config/...` | Config path relative to Axi_Schwarzschild | Use `../config/<name>.yaml` |
| `MGE file ... not found` | Wrong `template_directory` | Check templates dir exists at that path |
| `Directory not empty (rmtree)` | Cluster filesystem race; iterator crashed mid-cleanup | Just resubmit — directory is now clean |
| Failed `run.py` with no queue.txt | `clear_existing_models: true` + stale datfil/ | Set `clear_existing_models: false` or rm the dir manually, then resubmit |

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
fig.savefig('results/<name>/chi2_landscape.png', dpi=150, bbox_inches='tight')
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
| Check status | `python scripts/check_run.py <model_dir>` |
| Create config | `cp config/base.yaml config/new.yaml` + edit fields (⚠️ confirm with user) |
| Submit run | `run.py ../config/<name>.yaml` |
| Check jobs | `squeue -u wanght245001` |
| MPI verification | `srun -p test -n 2 --time=2 bash -c "source /etc/profile.d/module-profile.sh && module load openmpi/4.1.8 && which mpirun"` |
| Independent workers | `for i in $(seq 0 15); do sbatch -p test --mem 8G --wrap "source ... && conda activate schw && python scripts/schw_proc.py -p $i -q <queue>"; done` |
| Analyze results | `python scripts/analyze_results.py <model_dir> -o results/<name>/` |
| Output location | `results/<name>/` (default: same as model_dir) |
| Load grid | `Iterator(home_dir).get_model_list()` → `['dir','par','chi2']` |
| Best model dir | `model_list['dir'][0]` |
| Open model | `Model(dir, aperture_exprs=['_o','_s'])` |
| chi2 GP contour | `iter0.plot_chi2_grid()` → `(fig, axes)` |
| Kinematics | `Model.plot_kinematics()` |
| 1D slit | `Model.plot_kin1d(angle=90)` |
| 2D ellipsoids | `Model.plot_2D_anisotropy()` |
| Anisotropy | `Model.plot_anisotropy()` |
| Rerun NNLS | `rerun_nnls.py -i 0.02 -n 8 model_dir/*/` |
