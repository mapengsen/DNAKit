#!/usr/bin/env bash
# Update the self-hosted DNAKit web demo from this checkout:
#
#   packaging/huggingface_space/deploy.sh root@SERVER
#
# Copies src/dnakit/gui into the server's installed dnakit package (so the web page always matches
# the desktop app in this checkout), uploads the service files, and restarts the demo. The
# server must already be set up as described in packaging/README.md. One SSH connection is
# shared, so a password is asked once.
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 user@host" >&2
    exit 2
fi
HOST=$1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
CONTROL=$(mktemp -u "${TMPDIR:-/tmp}/dnakit-deploy.XXXXXX")
SSH=(ssh -o ControlMaster=auto -o "ControlPath=$CONTROL" -o ControlPersist=120 "$HOST")
trap 'ssh -o "ControlPath=$CONTROL" -O exit "$HOST" 2>/dev/null || true' EXIT

PKG=$("${SSH[@]}" '/opt/dnakit-web/venv/bin/python -c "import dnakit, pathlib; print(pathlib.Path(dnakit.__file__).parent)"')
echo "Updating dnakit.gui in $PKG"
tar -C "$ROOT/src/dnakit" --exclude=__pycache__ -czf - gui \
    | "${SSH[@]}" "rm -rf '$PKG/gui' && tar -xzf - -C '$PKG' && chmod -R a+rX '$PKG/gui'"

echo "Updating service files"
tar -C "$HERE" -czf - app.py requirements.txt README.md dnakit-web.service nginx-dnakit.conf \
    | "${SSH[@]}" 'tar -xzf - -C /opt/dnakit-web/app \
        && chmod -R a+rX /opt/dnakit-web/app \
        && cp /opt/dnakit-web/app/dnakit-web.service /etc/systemd/system/dnakit-web.service \
        && cp /opt/dnakit-web/app/nginx-dnakit.conf /etc/nginx/sites-available/dnakit \
        && systemctl daemon-reload && nginx -t -q && systemctl reload nginx \
        && systemctl restart dnakit-web'

echo "Waiting for the demo to start"
"${SSH[@]}" 'for i in $(seq 1 60); do
        curl -fs -o /dev/null http://127.0.0.1/ && { echo "OK: demo is serving"; exit 0; }
        sleep 2
    done
    journalctl -u dnakit-web -n 30 --no-pager; exit 1'
