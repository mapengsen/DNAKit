#!/usr/bin/env bash
# Update the self-hosted DNAKit web demo from this checkout:
#
#   packaging/huggingface_space/deploy.sh [--full] root@SERVER
#
# Copies src/dnakit/gui into the server's installed dnakit package (so the web page always matches
# the desktop app in this checkout), uploads the service files, and restarts the demo. With
# --full the whole src/dnakit package of this checkout replaces the installed one (the first
# time, the installed package is kept as dnakit.pypi-backup next to it for rollback). The
# server must already be set up as described in packaging/README.md. One SSH connection is
# shared, so a password is asked once.
set -euo pipefail

FULL=0
if [[ ${1:-} == --full ]]; then
    FULL=1
    shift
fi
if [[ $# -ne 1 ]]; then
    echo "usage: $0 [--full] user@host" >&2
    exit 2
fi
HOST=$1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
CONTROL=$(mktemp -u "${TMPDIR:-/tmp}/dnakit-deploy.XXXXXX")
SSH=(ssh -o ControlMaster=auto -o "ControlPath=$CONTROL" -o ControlPersist=120 "$HOST")
trap 'ssh -o "ControlPath=$CONTROL" -O exit "$HOST" 2>/dev/null || true' EXIT

PKG=$("${SSH[@]}" '/opt/dnakit-web/venv/bin/python -c "import dnakit, pathlib; print(pathlib.Path(dnakit.__file__).parent)"')
if [[ $FULL == 1 ]]; then
    echo "Replacing the dnakit package in $PKG with this checkout"
    # Stage the package as a wheel would install it: setup.py bundles the cataloged Skills into
    # dnakit/skills/_data at build time, so do the same (with the same checksum check) here.
    STAGE=$(mktemp -d "${TMPDIR:-/tmp}/dnakit-stage.XXXXXX")
    trap 'rm -rf "$STAGE"; ssh -o "ControlPath=$CONTROL" -O exit "$HOST" 2>/dev/null || true' EXIT
    tar -C "$ROOT/src" --exclude=__pycache__ -cf - dnakit | tar -xf - -C "$STAGE"
    python3 - "$ROOT" "$STAGE/dnakit/skills/_data" <<'PY'
import hashlib, json, shutil, sys
from pathlib import Path
root, target = Path(sys.argv[1]), Path(sys.argv[2])
catalog_path = root / "registry" / "catalog.json"
catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
shutil.rmtree(target, ignore_errors=True)
target.mkdir(parents=True)
(target / "catalog.json").write_bytes(catalog_path.read_bytes())
for skill in catalog["skills"]:
    for resource in skill["files"]:
        relative = Path(skill["id"]) / resource["path"]
        content = (root / "skills" / relative).read_bytes()
        if hashlib.sha256(content).hexdigest() != resource["sha256"]:
            sys.exit(f"Skill resource changed without rebuilding its catalog: {relative}")
        (target / relative).parent.mkdir(parents=True, exist_ok=True)
        (target / relative).write_bytes(content)
print(f"bundled {len(catalog['skills'])} Skills")
PY
    tar -C "$STAGE" -czf - dnakit \
        | "${SSH[@]}" "set -e; [ -e '$PKG.pypi-backup' ] || cp -a '$PKG' '$PKG.pypi-backup';
            rm -rf '$PKG.new' && mkdir '$PKG.new' && tar -xzf - -C '$PKG.new'
            rm -rf '$PKG' && mv '$PKG.new/dnakit' '$PKG' && rmdir '$PKG.new' && chmod -R a+rX '$PKG'"
else
    echo "Updating dnakit.gui in $PKG"
    tar -C "$ROOT/src/dnakit" --exclude=__pycache__ -czf - gui \
        | "${SSH[@]}" "rm -rf '$PKG/gui' && tar -xzf - -C '$PKG' && chmod -R a+rX '$PKG/gui'"
fi

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
