# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Achilles (A CHIcago Land Lepton Event Simulator) is a neutrino/electron-nucleus scattering event generator. It's a C++17/Fortran physics simulation with Python bindings via pybind11. The codebase uses CMake as its primary build system and scikit-build-core for pip-installable Python packaging.

## Build Commands

### CMake (C++ development)
```bash
cmake -B build -DACHILLES_ENABLE_TESTING=ON -DACHILLES_ENABLE_PYTHON=ON
cmake --build build -j$(nproc)
cmake --install build
```

### Python (pip install)
```bash
pip install .
```

### Key CMake options
- `-DACHILLES_ENABLE_TESTING=ON` — build test suite (Catch2 + Fortran tests)
- `-DACHILLES_ENABLE_PYTHON=ON` — build Python bindings
- `-DACHILLES_ENABLE_SHERPA=ON` — BSM calculations via Sherpa
- `-DACHILLES_ENABLE_ROOT=ON` — ROOT flux file support
- `-DACHILLES_ENABLE_HEPMC3=ON` — HepMC3 integration (default ON)
- `-DACHILLES_LOW_MEMORY=ON` — reduce memory usage during compilation

## Testing

```bash
# C++ tests (Catch2 v2.13.6 + trompeloeil for mocking)
./build/test/achilles-testsuite

# Run a single C++ test by name
./build/test/achilles-testsuite "[test-tag]"

# Python tests
pytest test/
```

## Formatting / Linting

- C++ formatting: clang-format v16.0.2 (config in `.clang-format`)
- Pre-commit hooks: `.pre-commit-config.yaml` (clang-format only)
- Run `pre-commit run --all-files` to format

## Architecture

### Language Boundaries
- **C++17** (`src/Achilles/`, `include/Achilles/`): Core physics engine — event generation, cascade simulation, integrators, nuclear models
- **Fortran** (`src/Achilles/fortran/`): Spectral functions, matrix element calculations (QE, RES models). Interfaces with C++ via explicit C bindings
- **Python** (`python/Achilles/`): Thin wrapper around the C++ `_achilles` extension module. Key exports: `generate_events()`, `initialize_logging()`, `__version__`

### Plugin System
Physics interactions, output formats, and cuts are loaded via a plugin architecture:
- `src/plugins/Cuts/` — event selection cuts
- `src/plugins/HepMC3/`, `src/plugins/NuHepMC/` — event output formats
- `src/plugins/Sherpa/` — BSM matrix elements
- `src/plugins/NucleonInteractions/` — cascade interaction models

### Configuration
Simulations are driven by YAML run cards (e.g., `run.yml`, `run_hnl.yml`). These configure:
- Beam/flux settings
- Nuclear targets and models
- Form factors (`FormFactors.yml`)
- Cascade parameters (`cascade.yml`)
- Output format and cuts

### Dependencies (fetched automatically via CPM)
fmt, spdlog, docopt, yaml-cpp, yaml-fortran, NuHepMC_CPPUtils, HighFive (HDF5), pybind11, Catch2

### Data Files
- `data/configurations/` — nuclear configuration files
- `data/default/` — default run card and form factor config
- `flux/` — neutrino flux configurations
- `densities/` — nuclear density files
- `UFO/` — dark neutrino portal BSM model definition

### Wheel data set (on-demand download)
The `data/` tree is ~227 MB — too large for a PyPI wheel/sdist (100 MB per-file limit),
so it is **not** bundled in the wheel. It is downloaded on demand into a cross-platform
user cache (via `platformdirs`) and located by both the Python module and the bundled
CLI through the `ACHILLES_DATA_DIR` env var (added to the C++ search path in
`Filesystem::AchillesPath`, `src/Achilles/System.cc`). Standalone CMake builds are
unaffected — `data/` still installs to `share/Achilles` (the wheel-only steps are
guarded by `SKBUILD`).

- `achilles-data status` / `get` / `path` — inspect and fetch the data set
  (`python/achilles/_data.py`). Downloads always prompt for consent (`-y` to skip).
  The bundled `achilles` CLI also offers to fetch on first run if data is missing.
- **Data is versioned independently of the code.** The required release tag is set
  once as `ACHILLES_DATA_RELEASE` in `CMakeLists.txt` and baked into the generated
  `python/achilles/_data_release.py` at configure time. Bump it only when `data/`
  changes.
- **To publish a new data set:** run `python scripts/package_data.py`, create a GitHub
  Release tagged `achilles-data-vN` (matching `ACHILLES_DATA_RELEASE`); the
  `achilles-data.yml` workflow attaches `achilles-data.tar.gz` + its `.sha256`.
