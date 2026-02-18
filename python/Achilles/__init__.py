from ._achilles import __version__, set_paths, generate_events
from importlib.resources import files
from pathlib import Path

root_dir = files("Achilles")
data_dir = root_dir / "data"
flux_dir = root_dir / "flux"

#TODO: Figure out a way for the plugin libs to be set or filled by the user using python?
libs_dir = root_dir / "lib"
libs_path = Path(libs_dir)
libs_path.mkdir(parents=True, exist_ok=True)

set_paths(root_dir, libs_dir, data_dir, flux_dir)

__all__ = ["__version__", "set_paths", "generate_events"]
