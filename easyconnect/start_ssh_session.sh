#!/bin/bash
# Establish SSH ControlMaster connection to cluster
# Subsequent SSH connections will reuse this session (no password prompts)

cd "$(dirname "$0")"
echo "=== Establishing SSH master connection to galaxy-login ==="
echo "Enter cluster password when prompted."
echo "Connection will persist for 4 hours (kept alive by ControlPersist)."
echo ""

ssh -o "ControlMaster=yes" -o "ControlPersist=4h" -N galaxy-login
