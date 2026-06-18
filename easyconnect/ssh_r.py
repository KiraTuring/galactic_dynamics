#!/usr/bin/env python3
"""Read-only SSH helper for cluster exploration.

Only allows safe/read commands (ls, du, cat, grep, stat, etc.).
Hard-blocks all write commands (rm, mv, mkdir, etc.) regardless of flags.

Usage:
  python3 ssh_r.py 'ls -la /some/dir'
  python3 ssh_r.py -t 600 'du -sh galaxy_models/'
"""

import sys, os

# Import core functions from ssh_helper
sys.path.insert(0, os.path.dirname(__file__))
from ssh_helper import _classify, ssh_run, _filter_output

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Read-only SSH helper for cluster inspection"
    )
    ap.add_argument("command", help="Command to run on cluster (read-only)")
    ap.add_argument("-t", "--timeout", type=int, default=300,
                    help="Command timeout in seconds (default: 300)")
    args = ap.parse_args()

    classification = _classify(args.command)

    if classification != 'safe':
        print(f"DENIED: ssh_r.py only allows read-only commands.\n"
              f"  Command: {args.command}\n"
              f"  Use 'python3 easyconnect/ssh_helper.py --exec ...' for write operations.",
              file=sys.stderr)
        sys.exit(2)

    output = ssh_run(args.command, timeout=args.timeout)
    _filter_output(output)
