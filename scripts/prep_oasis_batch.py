#!/usr/bin/env python3
"""Batch preprocess OASIS kinematics for a list of galaxies.

1. Decompress FITS.gz kinematics
2. Run create_kin_input() → aperture/bins/kinematics files
3. Create MGE table from sauron_mge.ecsv
"""
import os
import sys
import gzip
import shutil
import numpy as np
from astropy.io import ascii
from astropy import table as astropy_table

# Add data_prep dir to path
from pathlib import Path
_dp = Path(__file__).resolve().parent.parent / "Axi_Schwarzschild" / "data_prep"
sys.path.insert(0, str(_dp))

from generate_kin_input import create_kin_input
from data_combine import create_mge_table


DATA_RAW = os.path.abspath(os.path.join(os.path.dirname(__file__), "..",
                                         "data", "raw"))
DATA_PROC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..",
                                          "data", "processed"))
MGE_FILE = os.path.join(DATA_PROC, "sauron_mge.ecsv")


def decompress_kin(input_gz):
    """Decompress .fits.gz → .fits in the same directory."""
    output_fits = input_gz.replace(".fits.gz", ".fits")
    if os.path.exists(output_fits):
        print(f"  (already decompressed) {os.path.basename(output_fits)}")
        return output_fits
    with gzip.open(input_gz, "rb") as f_in:
        with open(output_fits, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"  decompressed → {os.path.basename(output_fits)}")
    return output_fits


def extract_mge(galaxy):
    """Extract MGE rows for a galaxy from sauron_mge.ecsv."""
    mge_all = ascii.read(MGE_FILE)
    mask = mge_all["galaxy"] == galaxy
    rows = mge_all[mask]
    if len(rows) == 0:
        raise ValueError(f"No MGE data found for {galaxy} in {MGE_FILE}")
    return rows["I"].tolist(), rows["sigma"].tolist(), rows["q"].tolist()


def prep_galaxy(galaxy):
    print(f"\n{'='*60}")
    print(f"Processing {galaxy}")
    print(f"{'='*60}")

    oasis_dir = os.path.join(DATA_RAW, galaxy, "oasis")
    out_dir = os.path.join(DATA_PROC, galaxy)
    os.makedirs(out_dir, exist_ok=True)

    # 1. Decompress kinematics
    print("\n[1/3] Decompressing kinematics...")
    kin_gz = os.path.join(oasis_dir, f"kinematics_oasis_{galaxy}.fits.gz")
    if not os.path.exists(kin_gz):
        print(f"  WARNING: {kin_gz} not found, skipping decompress.")
        kin_fits = None
    else:
        kin_fits = decompress_kin(kin_gz)

    # 2. Generate aperture, bins, kinematics
    print("\n[2/3] Generating aperture/bins/kinematics...")
    if kin_fits and os.path.exists(kin_fits):
        create_kin_input(
            galaxy=galaxy,
            file=kin_fits,
            dyn_model_dir=out_dir + os.sep,
            expr="_o",
            kin_input="OASIS",
            fit_PA=True,
            angle_deg=None,
            plot=True,
            min_gh_err=0.01,
        )
    else:
        print("  SKIPPED (no kinematics file)")

    # 3. Create MGE table
    print("\n[3/3] Creating MGE table...")
    lum, sigma_list, q_obs = extract_mge(galaxy)
    create_mge_table(
        dir=out_dir + os.sep,
        lum=lum,
        sigma=sigma_list,
        q_obs=q_obs,
    )
    print(f"  MGE: {len(lum)} components")

    print(f"\nDone: {out_dir}")


if __name__ == "__main__":
    galaxies = ["NGC4552", "NGC5846"]
    for g in galaxies:
        prep_galaxy(g)
    print(f"\n{'='*60}")
    print("All done.")
