#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
#
# install_nuclear_model.sh — Compile a nuclear-model plugin source file into a
# shared library and install it where Achilles discovers plugins.
#
# This is a thin wrapper around build_plugin.sh: it just links the plugin
# against achilles::nuclear_models (instead of the default achilles::cuts) so
# the NuclearModel base class and factory are available. Everything else --
# locating the installed Achilles via find_package(Achilles), generating the
# CMakeLists, building, and installing libAchillesPlugin_<name>.so into the
# plugin directory -- is handled by build_plugin.sh.
#
# Generate a starting source file with the companion script:
#   scripts/new_nuclear_model.sh
#
# Usage:
#   scripts/install_nuclear_model.sh [-n NAME] [-o OUTPUT_DIR] SOURCE.cc
#
#   SOURCE.cc   Nuclear-model plugin source file to compile (required).
#   -n NAME     Plugin name; the library is libAchillesPlugin_<name>.so.
#               Defaults to the source file's base name with a trailing
#               "NuclearModel" stripped (e.g. MySpectralNuclearModel.cc -> MySpectral).
#   -o DIR      Directory to install the .so into (default: Achilles's plugin dir).
#
set -euo pipefail

NAME=""
OUTPUT_DIR=""

while getopts ":n:o:h" opt; do
    case "$opt" in
        n) NAME="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :)  echo "Option -$OPTARG requires an argument" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

die() { printf '\033[1;31m[install_nuclear_model] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

SOURCE="${1:-}"
[[ -n "$SOURCE" ]]   || die "No source file given. Usage: $0 [-n NAME] [-o DIR] SOURCE.cc"
[[ -f "$SOURCE" ]]   || die "Source file '$SOURCE' does not exist."

# Derive a plugin name from the file name if one was not supplied:
# <Name>NuclearModel.cc -> <Name>, otherwise just strip the extension.
if [[ -z "$NAME" ]]; then
    base="$(basename "$SOURCE")"
    base="${base%.*}"
    NAME="${base%NuclearModel}"
fi
[[ "$NAME" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] \
    || die "Derived plugin name '$NAME' is not a valid C++ identifier; pass one with -n."

# Delegate to build_plugin.sh, linking the nuclear-model library so the plugin
# can see achilles::NuclearModel and the factory it registers into.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_args=(-n "$NAME" -s "$SOURCE" -l "achilles::nuclear_models")
[[ -n "$OUTPUT_DIR" ]] && build_args+=(-o "$OUTPUT_DIR")

exec "$SCRIPT_DIR/build_plugin.sh" "${build_args[@]}"
