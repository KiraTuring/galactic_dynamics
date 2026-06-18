#!/usr/bin/env python3
import os, sys, pty, select, signal, time

PASSWORD = sys.argv[1]
PROXY_CMD = "nc -X 5 -x localhost:1080 %h %p"

def ssh_with_password(password, host, cmd, timeout=30):
    pid, fd = pty.fork()
    if pid == 0:
        args = ["ssh", "-o", f"ProxyCommand={PROXY_CMD}",
                "-o", "StrictHostKeyChecking=accept-new", "-o", "ServerAliveInterval=10",
                "-o", "PreferredAuthentications=keyboard-interactive,password",
                host, cmd]
        os.execvp("ssh", args)
    
    output = b""
    sent = False
    t0 = time.time()
    try:
        while time.time() - t0 < timeout:
            r, _, _ = select.select([fd], [], [], 0.3)
            if r:
                try:
                    data = os.read(fd, 4096)
                except OSError:
                    break
                if not data:
                    break
                output += data
                s = data.decode(errors='replace').lower()
                if not sent and ("password:" in s):
                    os.write(fd, (password + "\n").encode())
                    sent = True
                    time.sleep(0.5)
    finally:
        os.close(fd)
    return output.decode(errors='replace')

pubkey = open(os.path.expanduser("~/.ssh/id_ed25519.pub")).read().strip()
cmd = f"mkdir -p ~/.ssh && chmod 700 ~/.ssh && echo '{pubkey}' >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys && echo KEY_OK"

print("=== Uploading SSH public key to cluster ===")
result = ssh_with_password(PASSWORD, "galaxy-login", cmd, timeout=25)
print(result)

if "KEY_OK" in result:
    print("\n=== Key uploaded! Testing key-based auth ===")
    r2 = ssh_with_password("wrong", "galaxy-login", "echo AUTH_OK", timeout=10)
    if "AUTH_OK" in r2:
        print("SUCCESS: Key-based authentication works!")
    else:
        print("Key auth test result:", r2)
else:
    print("Failed to detect success in output")
