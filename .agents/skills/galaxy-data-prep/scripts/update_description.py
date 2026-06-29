#!/usr/bin/env python3
"""
Update description.yaml for a galaxy. Designed to be both importable and CLI-callable.

Usage (CLI):
    python update_description.py NGC4552 --log "SAURON done" "PA diff=1.0 deg"
    python update_description.py NGC4552 --set data_quality "new text"
    python update_description.py NGC4552 \\
        --jam '{"name":"psf_free","chi2":2114,"q":0.695,"ratio":1.097,"lg_mbh":9.09,"lg_ml":0.594,"psf_source":"free"}'
    python update_description.py NGC4552 \\
        --jam '{"name":"default","chi2":2277}' --overwrite
    python update_description.py NGC4552 --jam-file psf_result.yaml

Usage (Python):
    from update_description import append_log, set_field, update_jam_model
    append_log("NGC4552", date="2026-01-01", step="PSF fit", notes="sigma=0.29")
    set_field("NGC4552", data_quality="...", redshift="0.003", overwrite=True)
    update_jam_model("NGC4552", name="default", chi2=2277, q=0.515, overwrite=True)
"""

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import yaml


def _resolve_root():
    """Resolve repo root by walking up from script dir until data/processed/ is found."""
    current = Path(__file__).resolve().parent
    for _ in range(10):
        if (current / "data" / "processed").is_dir():
            return current
        current = current.parent
    # Fallback: try cwd
    current = Path.cwd()
    for _ in range(10):
        if (current / "data" / "processed").is_dir():
            return current
        current = current.parent
    raise FileNotFoundError("Cannot find repo root (data/processed/ not found walking up)")


def _get_path(galaxy, root=None):
    if root is None:
        root = _resolve_root()
    return root / "data" / "processed" / galaxy / "description.yaml"


def _read(galaxy, root=None):
    p = _get_path(galaxy, root)
    if not p.exists():
        raise FileNotFoundError(f"{p} not found")
    with open(p) as f:
        return yaml.safe_load(f), p


def _write(data, galaxy, path):
    with open(path, "w") as f:
        f.write(f"# {galaxy} -- data description\n")
        yaml.dump(data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)


# ── public API ──────────────────────────────────────────────

def append_log(galaxy, *, step, notes, log_date=None, root=None):
    """Append a processing_log entry."""
    if log_date is None:
        log_date = date.today().isoformat()
    data, path = _read(galaxy, root)
    data.setdefault("processing_log", []).append(
        {"date": log_date, "step": step, "notes": notes}
    )
    _write(data, galaxy, path)
    print(f"  [{galaxy}] log appended: {step}")


def set_field(galaxy, root=None, overwrite=False, **kwargs):
    """Set top-level fields. e.g. set_field('NGC4552', data_quality='...', redshift='0.003', overwrite=True)"""
    if not kwargs:
        return
    data, path = _read(galaxy, root)
    if not overwrite:
        existing = [k for k in kwargs if k in data and data[k] is not None]
        if existing:
            print(f"  [{galaxy}] SKIPPED: fields already exist (use overwrite=True): {existing}")
            return
    data.update(kwargs)
    _write(data, galaxy, path)
    print(f"  [{galaxy}] updated: {', '.join(kwargs.keys())}")


def update_jam_model(galaxy, *, root=None, overwrite=False, **kwargs):
    """Add or replace a jam_models entry (matched by name). Requires overwrite=True to replace an existing entry."""
    name = kwargs.get("name")
    if not name:
        raise ValueError("Must specify 'name' for jam model entry")
    data, path = _read(galaxy, root)
    models = data.setdefault("jam_models", [])
    for i, m in enumerate(models):
        if m.get("name") == name:
            if not overwrite:
                print(f"  [{galaxy}] SKIPPED: model '{name}' already exists (use overwrite=True to replace)")
                return
            models[i] = kwargs
            _write(data, galaxy, path)
            print(f"  [{galaxy}] jam model replaced: {name}")
            return
    models.append(kwargs)
    _write(data, galaxy, path)
    print(f"  [{galaxy}] jam model appended: {name}")


# ── CLI ─────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Update description.yaml")
    parser.add_argument("galaxy", help="Galaxy name (e.g., NGC4552)")
    parser.add_argument("--log", nargs=2, metavar=("STEP", "NOTES"),
                        help="Append a processing_log entry")
    parser.add_argument("--date", help="Date for log entry (default: today)")
    parser.add_argument("--set", nargs=2, metavar=("KEY", "VALUE"),
                        action="append", help="Set a top-level field (repeatable)")
    parser.add_argument("--jam", metavar="JSON",
                        help="Add or replace jam model (append if new, replace if name exists + --overwrite)")
    parser.add_argument("--jam-file", metavar="YAML_PATH",
                        help="Read jam model params from a YAML file and append")
    parser.add_argument("--overwrite", action="store_true",
                        help="Allow overwriting existing jam model entries or description fields")
    args = parser.parse_args()

    did_something = False

    if args.log:
        log_date = args.date or date.today().isoformat()
        append_log(args.galaxy, step=args.log[0], notes=args.log[1], log_date=log_date)
        did_something = True

    if args.set:
        kv = {k: v for k, v in args.set}
        set_field(args.galaxy, overwrite=args.overwrite, **kv)
        did_something = True

    if args.jam:
        try:
            d = json.loads(args.jam)
        except json.JSONDecodeError as e:
            print(f"Error parsing --jam JSON: {e}")
            sys.exit(1)
        update_jam_model(args.galaxy, overwrite=args.overwrite, **d)
        did_something = True

    if args.jam_file:
        with open(args.jam_file) as f:
            d = yaml.safe_load(f)
        update_jam_model(args.galaxy, overwrite=args.overwrite, **d)
        did_something = True

    if not did_something:
        parser.print_help()


if __name__ == "__main__":
    main()
