#!/usr/bin/env python3
"""
Query NASA/IPAC Extragalactic Database (NED) for galaxy parameters:
redshift, morphology, coordinates, velocity, etc.

Designed to be both importable and CLI-callable.
Makes a single `query_object` call per galaxy (no secondary table queries).

Usage (CLI):
    python query_ned.py NGC3379
    python query_ned.py NGC4564 NGC5845
    python query_ned.py NGC3379 --yaml
    python query_ned.py NGC4564 --fill-description

Usage (Python):
    from query_ned import query_galaxy
    info = query_galaxy("NGC3379")
    print(info["redshift"])
"""

import argparse
import sys
import warnings
from pathlib import Path

import yaml


def _resolve_root():
    """Resolve repo root by walking up from script dir until data/processed/ is found."""
    current = Path(__file__).resolve().parent
    for _ in range(10):
        if (current / "data" / "processed").is_dir():
            return current
        current = current.parent
    current = Path.cwd()
    for _ in range(10):
        if (current / "data" / "processed").is_dir():
            return current
        current = current.parent
    return None


def query_galaxy(galaxy: str) -> dict:
    """
    Query NED for basic galaxy parameters.

    Returns a dict with keys:
        galaxy, ra_deg, dec_deg,
        redshift, redshift_flag,
        velocity_km_s,
        morphology, magnitude_filter,
        ned_url
    """
    result = {
        "galaxy": galaxy,
        "ra_deg": None, "dec_deg": None,
        "redshift": None, "redshift_flag": None,
        "velocity_km_s": None,
        "morphology": None,
        "magnitude_filter": None,
        "ned_url": f"https://ned.ipac.caltech.edu/byname?objname={galaxy}",
    }

    try:
        from astroquery.ipac.ned import Ned
    except ImportError:
        print(f"  WARNING: astroquery not available. pip install astroquery", file=sys.stderr)
        return result

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obj = Ned.query_object(galaxy)
    except Exception as exc:
        print(f"  [{galaxy}] NED query failed: {exc}", file=sys.stderr)
        return result

    if obj is None or len(obj) == 0:
        print(f"  [{galaxy}] NED: no results", file=sys.stderr)
        return result

    o = obj[0]
    try:
        result["ra_deg"] = float(o["RA"])
    except (KeyError, TypeError, ValueError):
        pass
    try:
        result["dec_deg"] = float(o["DEC"])
    except (KeyError, TypeError, ValueError):
        pass
    try:
        result["velocity_km_s"] = float(o["Velocity"])
    except (KeyError, TypeError, ValueError):
        pass
    try:
        result["morphology"] = str(o["Type"]).strip()
    except (KeyError, TypeError):
        pass
    try:
        result["redshift_flag"] = str(o["Redshift Flag"]).strip()
    except (KeyError, TypeError):
        pass
    try:
        result["magnitude_filter"] = str(o["Magnitude and Filter"]).strip()
    except (KeyError, TypeError):
        pass
    try:
        result["redshift"] = float(o["Redshift"])
    except (KeyError, TypeError, ValueError):
        pass

    return result


def format_table(galaxies: list[str]) -> str:
    """Pretty-print a table of results for multiple galaxies."""
    rows = [query_galaxy(g) for g in galaxies]
    header = f"{'Galaxy':<12s} {'z':>12s} {'flag':>8s} {'V':>8s} {'Morphology':<18s}"
    lines = [header, "-" * len(header)]
    for r in rows:
        z = f"{r['redshift']:.6f}" if r['redshift'] else "        --"
        zf = r['redshift_flag'] or "--"
        v = f"{r['velocity_km_s']:>7.0f}" if r['velocity_km_s'] else "    --"
        mor = r['morphology'] or "--"
        lines.append(f"{r['galaxy']:<12s} {z:>12s} {zf:>8s} {v:>8s} {mor:<18s}")
    return "\n".join(lines)


def to_yaml_results(galaxies: list[str]) -> str:
    """Return YAML string of results."""
    rows = [query_galaxy(g) for g in galaxies]
    if len(rows) == 1:
        return yaml.dump(rows[0], default_flow_style=False, sort_keys=False, allow_unicode=True)
    return yaml.dump(rows, default_flow_style=False, sort_keys=False, allow_unicode=True)


def fill_description(galaxy: str, root=None):
    """
    Update description.yaml for a galaxy with redshift from NED.

    Only writes fields currently marked 'TODO' or absent.
    """
    rv = query_galaxy(galaxy)
    if root is None:
        root = _resolve_root()
    if root is None:
        print(f"  ERROR: Cannot find repo root", file=sys.stderr)
        return

    yaml_path = root / "data" / "processed" / galaxy / "description.yaml"
    if not yaml_path.exists():
        print(f"  [{galaxy}] SKIP: {yaml_path} not found (run preprocessing first)")
        return

    with open(yaml_path) as f:
        data = yaml.safe_load(f) or {}

    changed = False
    if rv["redshift"] and data.get("redshift") in (None, "TODO"):
        data["redshift"] = f"{rv['redshift']:.6f}"
        changed = True
        print(f"  [{galaxy}] redshift = {data['redshift']}")

    if rv["morphology"] and data.get("type") in (None, "TODO"):
        data["type"] = rv["morphology"]
        changed = True
        print(f"  [{galaxy}] type = {data['type']}")

    if changed:
        with open(yaml_path, "w") as f:
            data["_updated_from_ned"] = True
            yaml.dump(data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
        print(f"  Updated: {yaml_path}")
    else:
        print(f"  [{galaxy}] nothing to update")


def main():
    parser = argparse.ArgumentParser(
        description="Query NED for galaxy parameters (redshift, morphology, etc.)")
    parser.add_argument("galaxies", nargs="+", help="Galaxy names (e.g. NGC3379 NGC4564)")
    parser.add_argument("--yaml", action="store_true", help="Output as YAML")
    parser.add_argument("--fill-description", action="store_true",
                        help="Update description.yaml for each galaxy with redshift & type from NED")
    args = parser.parse_args()

    if args.fill_description:
        for g in args.galaxies:
            fill_description(g)
    elif args.yaml:
        print(to_yaml_results(args.galaxies))
    else:
        print(format_table(args.galaxies))


if __name__ == "__main__":
    main()
