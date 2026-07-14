#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
#
# build_plugin.sh — Compile an Achilles plugin against an installed Achilles and
# drop the resulting shared library where Achilles looks for plugins.
#
# Achilles discovers plugins at runtime (src/plugins/Manager/PluginManager.cc):
# it scans every directory returned by GetPluginPaths() -- i.e. the entries in
# $ACHILLES_PLUGIN_PATH plus the install "lib" directory -- for files matching
# the regex  libAchillesPlugin.*\.so  and dlopen()s them. Each plugin must
# export an ExpectedVersion symbol (via the EXPORT_ACHILLES_VERSION() macro) and
# self-registers its component through static initialisation.
#
# This script:
#   1. Finds the Achilles source tree and its installed CMake package (the
#      AchillesConfig.cmake / AchillesTargets.cmake written by `make install`).
#   2. Writes a CMakeLists.txt that consumes Achilles via find_package(Achilles)
#      and links the imported achilles:: targets.
#   3. Configures + builds the plugin with CMake.
#   4. Installs the library into the plugin discovery directory.
#
# Usage:
#   scripts/build_plugin.sh [-n NAME] [-s SOURCE.cc] [-l LINK] [-o OUTPUT_DIR] [-r ACHILLES_ROOT]
#
#   -n NAME     Plugin name (default: Example). Library is libAchillesPlugin_<name>.so
#   -s SOURCE   Plugin C++ source file. If omitted, an example cut plugin is generated.
#   -l LINK     Imported Achilles target(s) to link (default: achilles::cuts). Use a
#               ';'-separated list for several, e.g. "achilles::nuclear_models".
#   -o DIR      Directory to install the .so into (default: auto-detected plugin dir).
#   -r DIR      Achilles source root (default: auto-detected from this script's location).
#
set -euo pipefail

# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
NAME="Example"
SOURCE=""
LINK="achilles::cuts"
OUTPUT_DIR=""
ACHILLES_ROOT=""

while getopts ":n:s:l:o:r:h" opt; do
    case "$opt" in
        n) NAME="$OPTARG" ;;
        s) SOURCE="$OPTARG" ;;
        l) LINK="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) ACHILLES_ROOT="$OPTARG" ;;
        h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :)  echo "Option -$OPTARG requires an argument" >&2; exit 1 ;;
    esac
done

log() { printf '\033[1;34m[build_plugin]\033[0m %s\n' "$*" >&2; }
die() { printf '\033[1;31m[build_plugin] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

# --------------------------------------------------------------------------- #
# 1. Locate the Achilles source root (top-level CMakeLists.txt + include/Achilles)
# --------------------------------------------------------------------------- #
if [[ -z "$ACHILLES_ROOT" ]]; then
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    dir="$script_dir"
    while [[ "$dir" != "/" ]]; do
        if [[ -f "$dir/CMakeLists.txt" && -d "$dir/include/Achilles" ]]; then
            ACHILLES_ROOT="$dir"
            break
        fi
        dir="$(dirname "$dir")"
    done
fi
[[ -n "$ACHILLES_ROOT" && -f "$ACHILLES_ROOT/CMakeLists.txt" ]] \
    || die "Could not locate the Achilles source root. Pass it with -r."
log "Achilles source root:      $ACHILLES_ROOT"

# --------------------------------------------------------------------------- #
# 2. Find the installed Achilles CMake package. `make install` writes
#    <prefix>/lib/cmake/Achilles/{AchillesConfig,AchillesTargets}.cmake; we pass
#    that directory to find_package(Achilles) via Achilles_DIR, and the install
#    prefix on CMAKE_PREFIX_PATH so the config's find_dependency(fmt/spdlog)
#    calls resolve against the bundled fmt/spdlog packages.
# --------------------------------------------------------------------------- #
ACHILLES_CONFIG="$(find "$ACHILLES_ROOT" -name 'AchillesConfig.cmake' \
                   -path '*cmake/Achilles*' 2>/dev/null | head -n1 || true)"
[[ -n "$ACHILLES_CONFIG" ]] || die \
    "Could not find an installed AchillesConfig.cmake under $ACHILLES_ROOT. Run 'make install' first."

ACHILLES_DIR="$(cd "$(dirname "$ACHILLES_CONFIG")" && pwd)"
# <prefix>/lib/cmake/Achilles  ->  <prefix>
PREFIX="$(cd "$ACHILLES_DIR/../../.." && pwd)"
log "Achilles CMake package:    $ACHILLES_DIR"
log "Achilles install prefix:   $PREFIX"

[[ -f "$ACHILLES_DIR/AchillesTargets.cmake" ]] || die \
    "AchillesConfig.cmake found but AchillesTargets.cmake is missing next to it -- the install did not export its targets."

# --------------------------------------------------------------------------- #
# 3. Determine where to install the finished plugin (plugin discovery dir).
#    Achilles scans <install-prefix>/lib (PathVariables::LibsPath) by default.
# --------------------------------------------------------------------------- #
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$PREFIX/lib"
fi
mkdir -p "$OUTPUT_DIR"
log "Plugin install directory:  $OUTPUT_DIR"

# --------------------------------------------------------------------------- #
# 4. Prepare the plugin source (generate an example if none supplied)
# --------------------------------------------------------------------------- #
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/achilles-plugin-${NAME}.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

if [[ -z "$SOURCE" ]]; then
    SOURCE="$WORK_DIR/${NAME}Plugin.cc"
    log "No source given; generating an example single-particle cut plugin."
    cat > "$SOURCE" <<EOF
// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Auto-generated example Achilles plugin: a single-particle cut named "${NAME}".
// Registers itself with the cut factory via static initialisation and exports
// the Achilles version the PluginManager checks before loading.
#include "Achilles/OneParticleCuts.hh"
#include "Achilles/Version.hh"

EXPORT_ACHILLES_VERSION()

namespace achilles {

ONE_PARTICLE_CUT(${NAME});

bool ${NAME}Cut::MakeCut(const FourVector &) const {
    return true;
}

} // namespace achilles
EOF
else
    [[ -f "$SOURCE" ]] || die "Plugin source '$SOURCE' does not exist."
    SOURCE="$(cd "$(dirname "$SOURCE")" && pwd)/$(basename "$SOURCE")"
fi
log "Plugin source:             $SOURCE"

# --------------------------------------------------------------------------- #
# 5. Write the CMakeLists.txt for the plugin.
#    find_package(Achilles) imports achilles::cuts (and its transitive deps:
#    achilles::utilities, the public headers, spdlog, yaml-cpp, ...), so the
#    plugin needs to name only the components it uses.
# --------------------------------------------------------------------------- #
CMAKE_FILE="$WORK_DIR/CMakeLists.txt"
cat > "$CMAKE_FILE" <<EOF
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Auto-generated by scripts/build_plugin.sh -- builds a single Achilles plugin
# against an installed Achilles via find_package(Achilles).
cmake_minimum_required(VERSION 3.16)
project(AchillesPlugin_${NAME} LANGUAGES CXX)

find_package(Achilles REQUIRED)

# Achilles loads any file named libAchillesPlugin*.so from its plugin path, so
# the target name must carry that prefix. MODULE == dlopen-only shared object.
add_library(AchillesPlugin_${NAME} MODULE "${SOURCE}")

# The linked achilles:: target pulls in the public headers, C++ standard and the
# rest of the library API transitively (e.g. achilles::cuts for a cut plugin,
# achilles::nuclear_models for a nuclear-model plugin).
target_link_libraries(AchillesPlugin_${NAME} PRIVATE ${LINK})

# Emit exactly libAchillesPlugin_${NAME}.so (no extra target-name mangling).
set_target_properties(AchillesPlugin_${NAME} PROPERTIES
    PREFIX "lib"
    OUTPUT_NAME "AchillesPlugin_${NAME}")
EOF
log "Wrote CMakeLists.txt:      $CMAKE_FILE"

# --------------------------------------------------------------------------- #
# 6. Configure + build
# --------------------------------------------------------------------------- #
CMAKE_BUILD="$WORK_DIR/build"
log "Configuring (find_package Achilles)..."
cmake -S "$WORK_DIR" -B "$CMAKE_BUILD" \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DAchilles_DIR="$ACHILLES_DIR" \
      -DCMAKE_PREFIX_PATH="$PREFIX" >&2
log "Building..."
cmake --build "$CMAKE_BUILD" --parallel >&2

ARTIFACT="$(find "$CMAKE_BUILD" -name "libAchillesPlugin_${NAME}.so" | head -n1)"
[[ -n "$ARTIFACT" ]] || die "Build succeeded but the plugin library was not produced."

# --------------------------------------------------------------------------- #
# 7. Install into the plugin discovery directory
# --------------------------------------------------------------------------- #
install -m 0755 "$ARTIFACT" "$OUTPUT_DIR/"
DEST="$OUTPUT_DIR/libAchillesPlugin_${NAME}.so"
log "Installed plugin ->        $DEST"
echo "$DEST"
