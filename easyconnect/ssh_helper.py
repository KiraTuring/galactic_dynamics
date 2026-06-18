#!/usr/bin/env python3
"""SSH helper for password-only cluster access via SOCKS5 proxy.

Usage:
  python3 ssh_helper.py <command>            safe commands auto-run; write commands BLOCKED
  python3 ssh_helper.py --exec <command>     allow write/destructive commands
  python3 ssh_helper.py -t 600 <command>     custom timeout (default: 300s)
  python3 ssh_helper.py --push <local> <remote>   upload to cluster
  python3 ssh_helper.py --pull <remote> <local>   download from cluster
  python3 ssh_helper.py --download <remote> <local>  pull single file

Permission model:
  Read commands (ls, du, cat, stat, etc.)  -> auto
  Write commands (rm, mv, mkdir, etc.)     -> require --exec
  No --exec with write command             -> blocked with error
"""

import os, sys, subprocess, time, pty, select, argparse, termios, tempfile

PASSWORD = os.environ.get("SSH_PASSWORD") or "R3w44CWWc*GATbk"
HOST = "galaxy-login"

# --------------- command classification ---------------

WRITE_PREFIXES = {
    'rm', 'mv', 'cp', 'mkdir', 'rmdir', 'touch', 'chmod', 'chown', 'ln',
    'tar', 'gzip', 'gunzip', 'zip', 'unzip', 'dd', 'shred', 'truncate',
    'mount', 'umount', 'kill', 'pkill', 'killall', 'reboot', 'shutdown',
    'sbatch', 'scancel', 'module', 'pip', 'pip3', 'rsync', 'scp',
    'systemctl', 'service', 'fdisk', 'mkfs', 'tee',
    'python', 'python3', 'conda',
}

def _classify(cmd):
    """Classify 'safe' or 'write'. Checks ALL chained commands (&& || ; |)."""
    # Split by all shell chain operators into independent commands
    segments = [cmd]
    for sep in ('&&', '||', ';', '|'):
        expanded = []
        for s in segments:
            expanded.extend(s.split(sep))
        segments = expanded

    # If any segment starts with a write command, the whole command is write
    for seg in segments:
        tokens = seg.strip().split()
        if not tokens:
            continue
        first = tokens[0].split('/')[-1]
        if first in WRITE_PREFIXES:
            return 'write'

    # Check redirects (write to files)
    if '>' in cmd:
        return 'write'

    return 'safe'


# --------------- ControlMaster management ---------------

# ssh_config already has ControlMaster auto for galaxy-login.
# But pty.fork() prevents the master socket from persisting.
# We use SSH_ASKPASS trick to start a background master (-MNf)
# which survives the Python process and enables fast subprocess reuse.


def _ensure_master():
    """Start a background ControlMaster via SSH_ASKPASS if not already alive."""
    if _master_alive():
        return True

    # Write password to a temporary script
    fd, script_path = tempfile.mkstemp(prefix="ssh_askpass_", suffix=".sh")
    with os.fdopen(fd, 'w') as f:
        f.write("#!/bin/bash\necho '%s'\n" % PASSWORD)
    os.chmod(script_path, 0o700)

    try:
        env = os.environ.copy()
        env["SSH_ASKPASS"] = script_path
        env["SSH_ASKPASS_REQUIRE"] = "force"
        subprocess.run(
            ["setsid", "ssh", "-MNf", "-o", "ControlPersist=4h", HOST],
            env=env, capture_output=True, timeout=15
        )
    finally:
        try:
            os.remove(script_path)
        except OSError:
            pass

    # Wait for socket to appear
    for _ in range(20):
        if _master_alive():
            return True
        time.sleep(0.2)
    return False


def _master_alive():
    """Check if master connection is responsive."""
    return subprocess.run(
        ["ssh", "-O", "check", HOST],
        capture_output=True, timeout=3
    ).returncode == 0


# --------------- command execution ---------------

def _pty_run(cmd, timeout=300):
    """Execute command via pty. First call auto-creates ControlMaster (via ssh_config)."""
    pid, fd = pty.fork()
    if pid == 0:
        os.execvp("ssh", ["ssh", HOST, cmd])

    # Disable echo to prevent password leak
    try:
        attr = termios.tcgetattr(fd)
        attr[3] &= ~termios.ECHO
        termios.tcsetattr(fd, termios.TCSANOW, attr)
    except termios.error:
        pass

    output = b""
    t0 = time.time()
    while time.time() - t0 < timeout:
        r, _, _ = select.select([fd], [], [], 0.3)
        if r:
            try:
                data = os.read(fd, 8192)
            except OSError:
                break
            if not data:
                break
            output += data
            if b"assword:" in output:
                os.write(fd, (PASSWORD + "\n").encode())

    os.close(fd)
    os.waitpid(pid, 0)
    return output.decode(errors='replace')


def ssh_run(cmd, timeout=300):
    """Execute command on cluster; prefers master reuse, falls back to pty."""
    # Master not alive? Try starting one (fast SSH_ASKPASS approach)
    if not _master_alive():
        _ensure_master()

    if _master_alive():
        result = subprocess.run(
            ["ssh", HOST, cmd],
            capture_output=True, text=True, timeout=timeout
        )
        if result.stderr:
            filtered = '\n'.join(l for l in result.stderr.split('\n')
                                 if 'Identity added' not in l and 'Agent pid' not in l
                                 and 'Pseudo-terminal' not in l and 'Connection to' not in l)
            if filtered.strip():
                sys.stderr.write(filtered.strip() + "\n")
        return result.stdout

    # Fallback: pty (will also create master, but slower)
    return _pty_run(cmd, timeout)


def _filter_output(text):
    """Remove noise lines from SSH output, redact password."""
    for line in text.split('\n'):
        skip = False
        for noise in ('Identity added', 'Agent pid', 'assword:',
                      'Pseudo-terminal', 'Connection to',
                      'ControlSocket', 'disabling multiplexing', PASSWORD):
            if noise in line:
                skip = True
                break
        if not skip and line.strip():
            print(line)


# --------------- CLI ---------------

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="SSH helper for HPC Westlake cluster. Write commands require --exec."
    )
    ap.add_argument("command", nargs="?", help="Command to run on cluster")
    ap.add_argument("-x", "--exec", dest="force_exec", action="store_true",
                    help="Allow write/destructive commands")
    ap.add_argument("-t", "--timeout", type=int, default=300,
                    help="Command timeout in seconds (default: 300)")
    ap.add_argument("--push", nargs=2, metavar=("LOCAL", "REMOTE"),
                    help="rsync LOCAL to cluster:REMOTE")
    ap.add_argument("--pull", nargs=2, metavar=("REMOTE", "LOCAL"),
                    help="Pull from cluster:REMOTE to LOCAL")
    ap.add_argument("--download", nargs=2, metavar=("REMOTE", "LOCAL"),
                    help="Download a file from cluster via ControlMaster (noise-free)")
    args = ap.parse_args()

    # --- sync modes ---
    for flag, direction in [('push', 'push'), ('pull', 'pull')]:
        paths = getattr(args, flag, None)
        if not paths:
            continue
        if not _master_alive():
            _ensure_master()
        src, dst = paths
        ssh_cmd = "ssh -o ControlMaster=auto"
        if flag == 'pull':
            rsync_cmd = f"rsync -avz --links -e '{ssh_cmd}' {HOST}:{src} {dst}"
        else:
            rsync_cmd = f"rsync -avz --links -e '{ssh_cmd}' {src} {HOST}:{dst}"
        print(f"[rsync] {rsync_cmd}")
        subprocess.run(rsync_cmd, shell=True, check=True)
        sys.exit(0)

    # --- download mode ---
    if args.download:
        if not _master_alive():
            _ensure_master()
        remote, local = args.download
        result = subprocess.run(
            ["ssh", HOST, "cat", remote],
            capture_output=True, timeout=120
        )
        # Strip remote noise lines from beginning (Agent pid / Identity added)
        data = result.stdout
        while data:
            nl = data.find(b'\n')
            if nl < 0 or nl > 100:
                break
            line = data[:nl]
            if line.startswith(b'Agent pid') or line.startswith(b'Identity added'):
                data = data[nl+1:]
            else:
                break
        with open(local, 'wb') as f:
            f.write(data)
        print(f"Downloaded {len(data)} bytes -> {local}")
        sys.exit(0)

    if not args.command:
        ap.print_help()
        sys.exit(1)

    cmd = args.command
    classification = _classify(cmd)

    if classification == 'write' and not args.force_exec:
        print(f"BLOCKED: Destructive command requires --exec (-x).\n"
              f"  Command: {cmd}\n"
              f"  Rerun with: python3 easyconnect/ssh_helper.py --exec '{cmd}'",
              file=sys.stderr)
        sys.exit(2)

    output = ssh_run(cmd, timeout=args.timeout)
    _filter_output(output)
