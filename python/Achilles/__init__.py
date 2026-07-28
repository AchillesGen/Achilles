from ._achilles import __version__, set_paths, initialize_logging, generate_events, verbosity
from pathlib import Path

# Register the fortran models
_achilles.register_fortran_models()

# Resolve the package directory from __file__ rather than importlib.resources.as_file:
# for an editable install files("Achilles") is a MultiplexedPath and as_file() extracts
# it to a *temporary* directory that is deleted when the with-block exits, leaving
# set_paths() pointing at paths that no longer exist (the C++ plugin loader then fails to
# open Achilles/lib). __file__ is a real, persistent path in both editable and wheel installs.
root_dir = Path(__file__).resolve().parent
data_dir = root_dir / "data"
flux_dir = root_dir / "flux"

# For editable installs, data/ and flux/ live at the repo root
# (two levels up from python/Achilles/), not inside the package directory.
if not data_dir.is_dir():
    _repo_root = root_dir.parent.parent
    if (_repo_root / "data").is_dir():
        root_dir = _repo_root
        data_dir = root_dir / "data"
        flux_dir = root_dir / "flux"

#TODO: Figure out a way for the plugin libs to be set or filled by the user using python?
libs_dir = root_dir / "lib"
libs_dir.mkdir(parents=True, exist_ok=True)

set_paths(root_dir, libs_dir, data_dir, flux_dir)

__all__ = ["__version__", "set_paths", "initialize_logging", "generate_events", "verbosity"]
