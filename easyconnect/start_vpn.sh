#!/bin/bash
# Start EasyConnect VPN (VNC mode)
# Connect with VNC Viewer to localhost:5901 (password: opencode)

cd "$(dirname "$0")"
echo "=== EasyConnect VPN ==="
echo "Starting container..."
echo "Connect VNC Viewer to localhost:5901 (password: opencode)"
echo "After login, SOCKS5 proxy available at localhost:1080"
echo ""
docker compose up -d
echo ""
echo "To view logs:  docker compose logs -f"
echo "To stop:       docker compose down"
