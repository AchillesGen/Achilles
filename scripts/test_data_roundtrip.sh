#!/bin/bash
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# End-to-end check of the on-demand data-set pipeline, entirely offline.
#
# Reproduces what the `data-release.yml` workflow uploads and what a
# `pip install achilles` user downloads, without touching GitHub or `main`:
#
#   1. Build the release assets with scripts/package_data.py (the workflow's
#      "Build data assets" step).
#   2. Serve dist/ over a throwaway HTTP server, standing in for the release
#      asset host.
#   3. Run the real downloader (achilles._data) against it via ACHILLES_DATA_URL
#      / ACHILLES_DATA_DIR, and confirm it downloads, verifies the SHA-256, and
#      extracts the data/ tree.
#   4. Negative test: corrupt the served tarball and confirm the checksum guard
#      rejects it.
#
# Usage: scripts/test_data_roundtrip.sh

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Downloader command: prefer the installed console script (the real user path),
# but fall back to running _data.py directly. The direct path avoids importing
# the achilles package __init__, so it works even in a source tree where the
# compiled achilles._achilles extension has not been built.
if command -v achilles-data >/dev/null 2>&1; then
    data_cmd=(achilles-data)
else
    data_cmd=(python "$repo_root/python/achilles/_data.py")
fi

work_dir="$(mktemp -d)"
dist_dir="$work_dir/dist"
cache_dir="$work_dir/cache"
server_pid=""

cleanup() {
    [[ -n "$server_pid" ]] && kill "$server_pid" 2>/dev/null || true
    rm -rf "$work_dir"
}
trap cleanup EXIT

pass() { printf '\033[32mPASS\033[0m %s\n' "$1"; }
fail() { printf '\033[31mFAIL\033[0m %s\n' "$1"; exit 1; }

# --- 1. Build the release assets (workflow's build step) --------------------
echo "== Building data assets =="
python "$repo_root/scripts/package_data.py" --data-dir "$repo_root/data" --out-dir "$dist_dir"
[[ -f "$dist_dir/achilles-data.tar.gz" ]]        || fail "tarball not produced"
[[ -f "$dist_dir/achilles-data.tar.gz.sha256" ]] || fail "checksum sidecar not produced"
pass "assets built"

# --- 2. Serve dist/ like a release asset host -------------------------------
port="$(python -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1]); s.close()')"
( cd "$dist_dir" && exec python -m http.server "$port" --bind 127.0.0.1 ) >/dev/null 2>&1 &
server_pid=$!

# Wait for the server to accept connections.
for _ in $(seq 1 50); do
    if python -c "import socket,sys; s=socket.socket(); sys.exit(s.connect_ex(('127.0.0.1',$port)))" 2>/dev/null; then
        break
    fi
    sleep 0.1
done

base_url="http://127.0.0.1:$port"
export ACHILLES_DATA_URL="$base_url/achilles-data.tar.gz"
export ACHILLES_DATA_DIR="$cache_dir"

# --- 3. Positive path: download, verify, extract ----------------------------
echo "== Downloading via the real consumer =="
"${data_cmd[@]}" -y get || fail "download/verify/extract failed"
"${data_cmd[@]}" status | grep -q "present" || fail "status did not report 'present'"
[[ -d "$cache_dir/data" ]] || fail "data/ tree not extracted into cache"
pass "download + checksum + extract"

# --- 4. Negative path: checksum guard rejects a corrupt tarball -------------
echo "== Verifying the checksum guard =="
printf 'corrupt' >> "$dist_dir/achilles-data.tar.gz"   # keep the original sidecar
bad_cache="$work_dir/cache-bad"
if ACHILLES_DATA_DIR="$bad_cache" "${data_cmd[@]}" -y get 2>/dev/null; then
    fail "corrupt tarball was accepted"
fi
[[ ! -d "$bad_cache/data" ]] || fail "corrupt tarball left extracted data behind"
pass "corrupt tarball rejected"

echo
echo "All data-roundtrip checks passed."
