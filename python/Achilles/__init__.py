from ._achilles import __version__, set_paths, generate_events
from importlib.resources import files, as_file
from pathlib import Path

# Register the fortran models
_achilles.register_fortran_models()

with as_file(files("Achilles")) as root_dir:
    root_dir = Path(root_dir)
    data_dir = root_dir / "data"
    flux_dir = root_dir / "flux"
    
    #TODO: Figure out a way for the plugin libs to be set or filled by the user using python?
    libs_dir = root_dir / "lib"
    libs_dir.mkdir(parents=True, exist_ok=True)
    
    set_paths(root_dir, libs_dir, data_dir, flux_dir)

__all__ = ["__version__", "set_paths", "generate_events"]
