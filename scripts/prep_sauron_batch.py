#!/usr/bin/env python3
"""SAURON kinematics preprocessing for NGC4552, NGC5846, NGC5813."""
import os
import sys
import yaml
import numpy as np
from pathlib import Path

# Add data_prep dir to path
_self = Path(__file__).resolve()
_dp = _self.parent.parent / "Axi_Schwarzschild" / "data_prep"
sys.path.insert(0, str(_dp))

from generate_kin_input import create_kin_input
from data_combine import offset_datafile, add_psf_to_datafile, combine_kin_file

ROOT = Path(__file__).resolve().parent.parent

galaxies = [
    {
        "name": "NGC4552",
        "release": "r1",
    },
    {
        "name": "NGC5846",
        "release": "r5",
    },
    {
        "name": "NGC5813",
        "release": "r3",
    },
]

# All PSF values are auto-read from HST fitting results below (no hardcoded defaults).


def read_oasis_pa(proc_dir):
    """Read PA from aperture_o.dat"""
    with open(proc_dir / "aperture_o.dat") as f:
        lines = f.readlines()
    # Format: header, then xmin ymin, xmax ymax, PA, nx ny
    pa = float(lines[3].strip())
    return pa


for g in galaxies:
    name = g["name"]
    release = g["release"]

    raw_dir = ROOT / "data" / "raw" / name
    proc_dir = ROOT / "data" / "processed" / name
    out_dir = str(proc_dir) + "/"

    cube_file = str(raw_dir / "sauron" / f"MS_{name}_{release}_C2D.fits")
    kin_file = str(raw_dir / "sauron" / f"{name}_{release}_idl.fits.gz")

    print(f"\n{'=' * 60}")
    print(f"  Processing {name} (release={release})")
    print(f"  Cube: {cube_file}")
    print(f"  Kin:  {kin_file}")
    print(f"  Out:  {out_dir}")
    print(f"{'=' * 60}")

    # Check files exist
    for f in [cube_file, kin_file]:
        if not os.path.exists(f):
            print(f"  ERROR: Missing {f}")
            sys.exit(1)

    # Read OASIS PA
    pa_o = read_oasis_pa(proc_dir)
    angle_deg = 90.0 - pa_o
    print(f"\n  OASIS PA: {pa_o:.2f}, using angle_deg={angle_deg:.2f}")

    # Define mask function: outside 4 arcsec from center
    def mask_func(xbin, ybin, kins):
        rbin = np.sqrt(xbin ** 2 + ybin ** 2)
        return rbin > 4

    # Step 1: Create _s kinematics input files
    print("\n  [1/5] Creating SAURON kinematics (_s)...")
    create_kin_input(
        name,
        [cube_file, kin_file],
        out_dir,
        expr="_s",
        fit_PA=True,
        kin_input="ATLAS3D",
        ngh=6,
        angle_deg=angle_deg,
        mask_func=mask_func,
        min_gh_err=0,
        gh_sys_err=[0.0, 0.0, 0.0, 0.0, 0.1, 0.1],
    )

    # Step 2: Apply sigma offset (SAURON sigma needs 1.05x correction)
    print("\n  [2/5] Applying sigma offset (x1.05)...")
    offset_datafile(out_dir + "gauss_hermite_kins_s.ecsv", offset={"sigma": 1.05})

    # Step 3: Add PSF metadata to both OASIS and SAURON data files.
    # Both are auto-read from their respective HST PSF fitting results.
    print("\n  [3/5] Adding PSF to OASIS and SAURON data files...")
    for aperture, subdir, label in [("_o", "psf_fit", "OASIS"), ("_s", "psf_fit_sauron", "SAURON")]:
        psf_file = ROOT / "results" / "JAM" / name / subdir / "psf_result.yaml"
        if not psf_file.exists():
            print(f"  ERROR: Missing {label} PSF result: {psf_file}")
            print(f"  Run: python JAM/scripts/fit_psf.py {name} --kin-input {'ATLAS3D' if aperture == '_s' else 'OASIS'}")
            sys.exit(1)
        psf = yaml.safe_load(open(psf_file))
        sigma = [psf["sigma1"], psf["sigma2"]]
        weight = [psf["weight1"], 1.0 - psf["weight1"]]
        print(f"  {label} PSF: sigma={sigma}, weight={weight}")
        add_psf_to_datafile(sigma, weight, out_dir + f"gauss_hermite_kins{aperture}.ecsv")

    # Step 4: Combine OASIS + SAURON
    print("\n  [4/5] Combining OASIS + SAURON kinematics...")
    combine_kin_file(out_dir, ["_o", "_s"])

    print(f"\n  {name} done!")

print("\n" + "=" * 60)
print("  All galaxies processed.")
print("=" * 60)
