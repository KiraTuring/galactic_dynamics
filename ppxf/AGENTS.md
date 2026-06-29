# pPXF Workbench

Stellar kinematic extraction from OASIS IFU datacubes using pPXF + MILES templates.

## Entrypoints

| File | Role |
|---|---|
| `run_ppxf.py` | **CLI**: full pPXF pipeline (read cube → bin → log-rebin → pPXF fit → save) |
| `plot_kin_comparison.py` | **CLI**: kinematic map comparison (pPXF vs reference OASIS) |
| `plot_scatter_comparison.py` | **CLI**: bin-by-bin scatter comparison (pPXF vs reference) |
| `ppxf.ipynb` | Interactive: exploration, parameter tuning, visualization of pPXF results |
| `ppxfprep.py` | Library: reading, binning, fitting, bootstrap, saving (imported by CLI and notebooks) |
| `ppxf_diagnostics.py` | Library: mock spectra, MC bias simulation, diagnostic plots (imported by notebooks) |
| `ppxf__test.py` | Smoke test (`ppxf.ppxf`, `ppxf.ppxf_util`, `ppxf.sps_util`) |

## CLI Usage

### `run_ppxf.py` — pPXF kinematic extraction

```bash
python run_ppxf.py NGC4621 --redshift 0.001438 [options]

# Required:
#   galaxy              Galaxy name (e.g. NGC4621)
#   --redshift FLOAT    Galaxy redshift

# Optional:
#   --cube PATH         OASIS FITS cube (default: auto-detect from ../data/raw/{GALAXY}/oasis/)
#   --miles PATH        MILES library dir (default: ./miles_lib/MILES_library_v9.1_FITS)
#   --target-sn FLOAT   Voronoi target S/N (default: 60)
#   --fwhm-gal FLOAT    Instrumental FWHM in Å (default: 5.4)
#   --moments INT       GH moments: 4 or 6 (default: 6)
#   --bias FLOAT        pPXF bias (default: 0.5)
#   --degree INT        Additive polynomial degree (default: 4)
#   --bootstrap INT     Bootstrap iterations for errors (0=skip, default: 100)
#   --output-dir PATH   Output directory (default: output/{GALAXY})
#   --trial-moments     Moments to test on total spectrum (default: 4 6)
#   --ref-bins [PATH]   Use reference binning from OASIS kinematics FITS (default: auto-detect)
```

With `--ref-bins`, uses the `BINNUM` column from a reference OASIS kinematics FITS file
instead of Voronoi binning. Output goes to `output/{GALAXY}_refbins/`.
Auto-detects `../data/raw/{GALAXY}/oasis/kinematics_oasis_{GALAXY}.fits*`.

**Important**: when using `--ref-bins --moments 6`, the comparison scripts (`plot_kin_comparison.py`,
`plot_scatter_comparison.py`) default to `h4` and non-`_refbins` paths. You must
specify all paths explicitly:

```bash
python plot_kin_comparison.py NGC4564 \
    --ecsv output/NGC4564_refbins/gauss_hermite_kins_h6.ecsv \
    --ref ../data/raw/NGC4564/oasis/kinematics_oasis_NGC4564.fits.gz \
    --output output/NGC4564_refbins/kin_comparison_h6.png
```

Output files in `output/{GALAXY}/` (or `output/{GALAXY}_refbins/` with `--ref-bins`):
- `ppxf_kins_mileslib_h{MOMENTS}_b{BIAS}.npz` — raw results (kin_list, dkin_list, x_gen, y_gen)
- `ppxf_kins_{GALAXY}_h{MOMENTS}.fits` — SAURON-format FITS kinematics
- `gauss_hermite_kins_h{MOMENTS}.ecsv` — GH kinematics table
- `description.yaml` — experiment log (auto-updated by `run_ppxf.py`, see below)

### `plot_kin_comparison.py` — Kinematic map comparison

```bash
python plot_kin_comparison.py NGC5813 [options]

# Options:
#   --ecsv PATH         Our pPXF ECSV file (default: output/{GALAXY}/gauss_hermite_kins_h4.ecsv)
#   --ref PATH          Reference OASIS kinematics FITS (default: auto-detect)
#   --cube PATH         OASIS cube for pixel coords (default: auto-detect)
#   --target-sn FLOAT   Voronoi target S/N for bin reconstruction (default: 40)
#   --systemic-vel FLOAT  Systemic velocity (default: median of our V)
#   --output PATH       Output PNG (default: output/{GALAXY}/kin_comparison.png)
```

Automatically detects h5/h6 columns → 6-column layout (V, σ, h3–h6).
Colorbar ranges symmetric about 0 for V and h3–h6; σ keeps positive range.

### `plot_scatter_comparison.py` — Bin-by-bin scatter comparison

```bash
python plot_scatter_comparison.py NGC5813 [options]

# Options:
#   --ecsv PATH         Our pPXF ECSV file (default: output/{GALAXY}_refbins/gauss_hermite_kins_h4.ecsv)
#   --ref PATH          Reference OASIS kinematics FITS (default: auto-detect)
#   --output PATH       Output PNG (default: output/{GALAXY}_refbins/scatter_comparison.png)
```

Automatically detects h5/h6 → 2×3 layout. V subtracted by median.
Dashed zero-lines for V and h3–h6 (not σ).

## Data Layout

```
miles_lib/         # MILES stellar library v9.1 (985 FITS spectra, from tar.gz)
output/{GALAXY}/   # Per-galaxy: .npz archives, FITS kinematics tables, ECSV kinematics
```

OASIS cubes are read from `../data/raw/{GALAXY}/oasis/` (auto-detected by `run_ppxf.py`).

## Conventions

- Galaxy names: uppercase NGC format (e.g., `NGC4621`)
- Kinematics columns: `v`, `dv`, `sigma`, `dsigma`, `h3`, `dh3`, `h4`, `dh4`, `h5`, `dh5`, `h6`, `dh6`, `is_good`, `n_gh`
- Instrumental FWHM for all OASIS data: `fwhm_gal = 5.4 Å`
- pPXF `bias=0.5` taken from McDermid et al. (2006)
- Never read FITS files into context — log paths only (per parent AGENTS.md)

### Known Galaxy Parameters

| Galaxy | redshift | N_bins (OASIS) | Notes |
|--------|----------|----------------|-------|
| NGC3379 | 0.003026 | 961 | Large — skip bootstrap (`--bootstrap 0`) |
| NGC4564 | 0.003809 | 370 | |
| NGC4621 | 0.001438 | 618 | |
| NGC5813 | 0.006540 | 459 | |
| NGC5845 | 0.005550 | 320 | |

Redshifts from NED (query via `.agents/skills/galaxy-data-prep/scripts/query_ned.py`).

## Key Dependencies

`ppxf`, `numpy`, `scipy`, `matplotlib`, `astropy`, `corner`, `vorbin`, `plotbin`

## Adding a New Galaxy

1. Run pPXF fitting (OASIS cube auto-detected from `../data/raw/{GALAXY}/oasis/`):
   ```bash
   python run_ppxf.py NGCXXXX --redshift Z --fwhm-gal FWHM --target-sn SN
   ```

## Notes

- CLI scripts are the primary interface for automated/agent-driven execution
- Notebooks are for interactive exploration, parameter tuning, and visualization only
- `ppxfprep.py` is the core library for the pipeline; `ppxf_diagnostics.py` contains mock/diagnostic tools for notebooks
- The parent project AGENTS.md (`../AGENTS.md`) covers overall project conventions and cluster access

## Performance

Bootstrap error estimation dominates runtime (100 iterations × N_bins × pPXF fit).

| Bins | Main fit | Bootstrap | Guideline |
|------|----------|-----------|-----------|
| ~320 | ~8s | ~420s | Default (`--bootstrap 100`) |
| ~370 | ~10s | ~600s | Default |
| ~460 | ~15s | ~800s | Borderline |
| ~620 | ~20s | ~1000s | Consider `--bootstrap 50` |
| ~960 | ~15s | ~1150s | **Must use `--bootstrap 0`** (or `setsid` for background) |

Bootstrap scales roughly O(N_bins) per iteration.

When running in background / redirecting to log files, use `python -u` to disable
stdout buffering — otherwise no progress output appears until the process exits.

## Debugging

- **Background processes with conda**: `source activate` may fail in subprocesses.
  Use the full Python path: `/path/to/miniconda3/envs/schw/bin/python3`.
  Use `setsid` to fully detach long-running processes — `nohup ... &` alone
  may be killed when the parent shell exits:

- **Reference FITS extension**: files are `.fits.gz` not `.fits`. Auto-detect handles
  this, but manually specified `--ref` must include the full filename with `.gz`.

- **Comparison script default paths**: `plot_kin_comparison.py` defaults to
  `output/{GALAXY}/gauss_hermite_kins_h4.ecsv`. When you used `--ref-bins` and/or
  `--moments 6`, the actual output is under `output/{GALAXY}_refbins/` with `h6`.
  Explicitly pass `--ecsv`, `--ref`, and `--output`.

## Experiment Logging

Each output directory contains a `description.yaml` with two sections:

- **`summary`** (manually maintained): free-text notes on findings, data provenance, and key results
- **`data`** (manually maintained): galaxy metadata (redshift, fwhm_gal, binning source, etc.)
- **`runs`** (auto-appended by `run_ppxf.py`): list of run records with timestamp, parameters, and output files

`run_ppxf.py` automatically appends a new entry to `runs` on each successful completion.
If `description.yaml` doesn't exist, it creates one with an empty `summary`/`data` and the new run record.
