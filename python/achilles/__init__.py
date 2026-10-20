# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

import os
from ._achilles import __version__, set_paths, initialize_logging, generate_events, verbosity
from importlib.resources import files, as_file
from pathlib import Path

from . import _data

# Register the fortran models
_achilles.register_fortran_models()

with as_file(files("achilles")) as root_dir:
    root_dir = Path(root_dir)
    data_dir = root_dir / "data"
    flux_dir = root_dir / "flux"

    # For editable installs, data/ and flux/ live at the repo root
    # (two levels up from python/achilles/), not inside the package directory.
    if not data_dir.is_dir():
        _repo_root = root_dir.parent.parent
        if (_repo_root / "data").is_dir():
            root_dir = _repo_root
            data_dir = root_dir / "data"
            flux_dir = root_dir / "flux"

    # Wheel installs ship no data/: fall back to the on-demand download cache
    # (populated by `achilles-data get`). data files are referenced with a "data/"
    # prefix, so the cache root is the parent of data_dir.
    if not data_dir.is_dir() and _data.data_root().is_dir():
        data_dir = _data.data_root()

    # Point the C++ search path (Filesystem::AchillesPath) at the cache so the module
    # and the bundled CLI resolve data identically. Set even when data isn't present
    # yet, so a later `achilles-data get` is picked up without re-importing.
    os.environ.setdefault("ACHILLES_DATA_DIR", str(_data.cache_dir()))

    #TODO: Figure out a way for the plugin libs to be set or filled by the user using python?
    libs_dir = root_dir / "lib"
    libs_dir.mkdir(parents=True, exist_ok=True)

    set_paths(root_dir, libs_dir, data_dir, flux_dir)

__all__ = ["__version__", "set_paths", "initialize_logging", "generate_events", "verbosity"]
