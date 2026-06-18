#!/usr/bin/env python3
"""
Generate a description.yaml stub for a galaxy's processed data directory.

Usage:
    python write_description.py NGC4552

Reads existing mge.ecsv, gauss_hermite_kins_o.ecsv, and aperture_o.dat from
data/processed/{galaxy}/ to fill in known fields. Unknown fields are marked
with YAML comments as TODO.

Embedded lookup table covers 28 galaxies from:
McDermid+2006, MNRAS 373, 906 (SAURON Project VIII).
"""

import os
import sys
import argparse
from datetime import date
from pathlib import Path

import numpy as np
import yaml
from astropy.io import ascii

# ---- McDermid+2006 Table 1 (28 OASIS galaxies) ----
_MCDERMID_TABLE = {
    "NGC1023": {"type": "SB0-(rs)", "fwhm": 1.41, "run": 2, "distMpc": 11.4},
    "NGC2549": {"type": "S0(r)sp", "fwhm": 0.91, "run": 1, "distMpc": 12.5},
    "NGC2695": {"type": "SAB0(s)", "fwhm": 0.91, "run": 2, "distMpc": 24.7},
    "NGC2699": {"type": "E:", "fwhm": 1.67, "run": 2, "distMpc": 24.7},
    "NGC2768": {"type": "E6:", "fwhm": 0.89, "run": 3, "distMpc": 21.4},
    "NGC2974": {"type": "E4", "fwhm": 0.99, "run": 2, "distMpc": 20.1},
    "NGC3032": {"type": "SAB0(r)", "fwhm": 0.92, "run": 2, "distMpc": 22.0},
    "NGC3379": {"type": "E1", "fwhm": 0.91, "run": 1, "distMpc": 10.3},
    "NGC3384": {"type": "SB0-(s):", "fwhm": 0.72, "run": 1, "distMpc": 11.1},
    "NGC3414": {"type": "S0 pec", "fwhm": 0.86, "run": 3, "distMpc": 25.2},
    "NGC3489": {"type": "SAB0+(rs)", "fwhm": 0.69, "run": 3, "distMpc": 12.0},
    "NGC3608": {"type": "E2", "fwhm": 1.13, "run": 1, "distMpc": 22.9},
    "NGC4150": {"type": "S0(r)?", "fwhm": 2.15, "run": 2, "distMpc": 13.7},
    "NGC4262": {"type": "SB0-(s)", "fwhm": 0.60, "run": 3, "distMpc": 17.2},
    "NGC4382": {"type": "S0+(s)pec", "fwhm": 1.11, "run": 2, "distMpc": 18.5},
    "NGC4459": {"type": "S0+(r)", "fwhm": 1.53, "run": 3, "distMpc": 16.4},
    "NGC4473": {"type": "E5", "fwhm": 0.80, "run": 3, "distMpc": 16.1},
    "NGC4486": {"type": "E0-1+pec", "fwhm": 0.57, "run": 1, "distMpc": 16.7},
    "NGC4526": {"type": "SAB0(s)", "fwhm": 0.94, "run": 1, "distMpc": 16.9},
    "NGC4552": {"type": "E0-1", "fwhm": 0.67, "run": 1, "distMpc": 15.9},
    "NGC4564": {"type": "E", "fwhm": 0.70, "run": 3, "distMpc": 17.8},
    "NGC4621": {"type": "E5", "fwhm": 0.86, "run": 3, "distMpc": 18.3},
    "NGC5198": {"type": "E1-2:", "fwhm": 0.84, "run": 3, "distMpc": 36.3},
    "NGC5308": {"type": "S0-", "fwhm": 1.91, "run": 5, "distMpc": 25.1},
    "NGC5813": {"type": "E1-2", "fwhm": 0.87, "run": 1, "distMpc": 32.2},
    "NGC5831": {"type": "E3", "fwhm": 0.95, "run": 3, "distMpc": 27.2},
    "NGC5845": {"type": "E:", "fwhm": 0.91, "run": 3, "distMpc": 25.0},
    "NGC5846": {"type": "E0-1", "fwhm": 0.58, "run": 3, "distMpc": 24.9},
    "NGC5982": {"type": "E3", "fwhm": 0.77, "run": 3, "distMpc": 41.9},
}


def _resolve_root():
    """Find project root relative to this script.
    Script: .agents/skills/galaxy-data-prep/scripts/write_description.py
    Root:   5 levels up.
    """
    return Path(__file__).resolve().parent.parent.parent.parent.parent


def _read_optional_ecsv(path):
    """Read an ECSV file, return None if it doesn't exist."""
    if not os.path.exists(path):
        return None
    return ascii.read(path)


def write_stub(galaxy, processed_dir=None, overwrite=False):
    """
    Write description.yaml stub for a galaxy.

    Parameters
    ----------
    galaxy : str
        Galaxy name (e.g. "NGC4552").
    processed_dir : str or Path, optional
        Path to data/processed/{galaxy}/. Defaults to inferred project root.
    overwrite : bool
        If False, skip if description.yaml already exists.
    """
    root = _resolve_root()
    if processed_dir is None:
        processed_dir = root / "data" / "processed" / galaxy
    else:
        processed_dir = Path(processed_dir)

    out_path = processed_dir / "description.yaml"
    if out_path.exists() and not overwrite:
        print(f"  (already exists, skipping) {out_path}")
        return

    today = date.today().isoformat()

    # ---- Build stub from McDermid table ----
    mc = _MCDERMID_TABLE.get(galaxy, {})
    stub = {
        "galaxy": galaxy,
        "type": mc.get("type", "TODO"),
        "distMpc": mc.get("distMpc", "TODO"),
        "redshift": "TODO",
        "fwhm_arcsec": mc.get("fwhm", "TODO"),
        "oasis_run": mc.get("run", "TODO"),
        "pixsize": 0.27,
    }

    # ---- MGE ----
    mge_tab = _read_optional_ecsv(processed_dir / "mge.ecsv")
    if mge_tab is not None:
        q = mge_tab["q"].data
        s = mge_tab["sigma"].data
        stub["mge"] = {
            "source": "sauron_mge.ecsv",
            "components": int(len(mge_tab)),
            "q_min": float(np.min(q)),
            "q_max": float(np.max(q)),
            "sigma_min_arcsec": float(np.min(s)),
            "sigma_max_arcsec": float(np.max(s)),
        }
    else:
        stub["mge"] = "TODO"

    # ---- Kinematics ----
    kin_tab = _read_optional_ecsv(processed_dir / "gauss_hermite_kins_o.ecsv")
    if kin_tab is not None:
        v = kin_tab["v"].data
        sig = kin_tab["sigma"].data
        stub["kinematics"] = {
            "bins": int(len(kin_tab)),
            "v_min_kms": float(np.min(v)),
            "v_max_kms": float(np.max(v)),
            "sigma_min_kms": float(np.min(sig)),
            "sigma_max_kms": float(np.max(sig)),
            "pa": "TODO",
        }
    else:
        stub["kinematics"] = "TODO"

    # Try to read PA from aperture file
    ap_file = processed_dir / "aperture_o.dat"
    if os.path.exists(ap_file):
        with open(ap_file) as f:
            lines = f.readlines()
            if len(lines) >= 4:
                raw_pa = float(lines[3].strip())
                fitted_pa = 90 - raw_pa
                if "kinematics" in stub and isinstance(stub["kinematics"], dict):
                    stub["kinematics"]["pa"] = round(fitted_pa, 1)

    # ---- Data quality (manual) ----
    if mc.get("fwhm"):
        fwhm_rank = sum(1 for v in _MCDERMID_TABLE.values()
                       if v.get("fwhm", 999) < mc["fwhm"]) + 1
        total = sum(1 for v in _MCDERMID_TABLE.values() if "fwhm" in v)
        stub["data_quality"] = (
            f"PSF FWHM {mc['fwhm']}\" (rank {fwhm_rank}/{total} in McDermid+2006). "
            "TODO: assess overall quality, note any issues."
        )
    else:
        stub["data_quality"] = "TODO: assess overall quality, note any issues."

    # ---- Processing log (empty) ----
    stub["processing_log"] = [{
        "date": today,
        "step": "stub generated",
        "notes": (
            f"mge: {len(mge_tab) if mge_tab is not None else 'N/A'} comps, "
            f"kin: {stub.get('kinematics', {}).get('bins', 'N/A') if isinstance(stub.get('kinematics'), dict) else 'N/A'} bins"
        ),
    }]

    # ---- JAM models (empty) ----
    stub["jam_models"] = []

    # ---- HST data ----
    raw_dir = root / "data" / "raw" / galaxy
    hst_files = []
    if raw_dir.exists():
        for f in sorted(raw_dir.glob("hst_*_drz.fits")):
            hst_files.append(f.name)
    stub["hst_images"] = hst_files if hst_files else "TODO"

    # ---- Write ----
    os.makedirs(processed_dir, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f"# {galaxy} — data description (auto-generated stub)\n")
        f.write(f"# Fields marked \"TODO\" require manual review.\n\n")
        yaml.dump(stub, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    print(f"  Wrote stub to {out_path}")


def append_log(galaxy, entry, processed_dir=None):
    """Append a log entry to an existing description.yaml. (delegates to update_description)"""
    from update_description import append_log as _append_log
    log = entry.get("processing_log", {})
    _append_log(
        galaxy,
        step=log.get("step", "unknown"),
        notes=log.get("notes", ""),
        log_date=log.get("date"),
        root=processed_dir
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate description.yaml stub for a galaxy"
    )
    parser.add_argument("galaxy", help="Galaxy name (e.g., NGC4552)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing description.yaml")
    args = parser.parse_args()

    write_stub(args.galaxy, overwrite=args.overwrite)
