#!/bin/bash
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# Build the Achilles wheel locally and smoke-test it in a fresh venv, without the
# PyPI upload cycle.
#
# The point is to reproduce a *relocated* wheel: `pip wheel` builds in a temp dir
# that is cleaned up before install, so any build-time absolute path baked into
# the binary (e.g. the scikit-build staging prefix) is gone at run time -- exactly
# the condition seen when installing from PyPI. Installing the built .whl into a
# fresh venv and launching the bundled `achilles` CLI exercises the same
# _cli.py -> compiled binary -> plugin scan path that fails on a bad install.
#
# The CMake build tree is cached under build/wheel, so after the first (slow)
# compile, iterating on a C++ change only recompiles what changed.
#
# Usage: scripts/build_and_test_wheel.sh [-- <args passed to the achilles CLI>]
#   e.g. scripts/build_and_test_wheel.sh -- --help

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# Arguments after `--` are forwarded to the achilles CLI smoke test. With no `--`,
# the CLI is launched with no arguments: its default "run mode" (runcard=run.yml),
# which is the path that exercises the relocated-wheel data-resolution fixes.
cli_args=()
if [[ "${1:-}" == "--" ]]; then
    shift
    cli_args=("$@")
fi

# Build and install in a throwaway dir so we never touch (or get confused by) any
# real wheels already sitting in wheelhouse/.
scratch="$(mktemp -d)"
wheel_out="$scratch/wheelhouse"
venv_dir="$scratch/wheeltest"

cleanup() { rm -rf "$scratch"; }
trap cleanup EXIT

pass() { printf '\033[32mPASS\033[0m %s\n' "$1"; }
fail() { printf '\033[31mFAIL\033[0m %s\n' "$1"; exit 1; }

# --- 1. Ensure the build backend is importable for --no-build-isolation -------
python -c "import scikit_build_core" 2>/dev/null || {
    echo "Installing scikit-build-core (needed for --no-build-isolation)"
    pip install "scikit-build-core>=0.5"
}

# --- 2. Build the wheel (cached CMake tree -> fast incremental rebuilds) -------
echo "== Building wheel =="
pip wheel . --no-deps --no-build-isolation \
    --config-settings=build-dir=build/wheel \
    -w "$wheel_out"

wheel="$(ls -1t "$wheel_out"/*.whl 2>/dev/null | head -n1)"
[[ -n "$wheel" ]] || fail "no wheel produced in $wheel_out"
pass "built $(basename "$wheel")"

# --- 3. Install the built artifact into a fresh, isolated venv -----------------
echo "== Installing into a fresh venv =="
python -m venv "$venv_dir"
"$venv_dir/bin/pip" install --quiet "$wheel"
pass "installed into $venv_dir"

# --- 4. Smoke-test the bundled CLI (the relocated-binary failure path) ---------
# Run from a scratch dir so a local ./data or ./lib can't mask a bad install.
echo "== Smoke-testing the achilles CLI: achilles ${cli_args[*]} =="
if ( cd "$(dirname "$venv_dir")" && "$venv_dir/bin/achilles" "${cli_args[@]}" ); then
    pass "CLI ran without terminating"
else
    fail "CLI exited non-zero (see output above)"
fi

echo
echo "Wheel build + install + CLI smoke test passed."
