"""Console-script launcher for the bundled native ``achilles`` executable.

The compiled C++ ``achilles`` binary is installed inside the package directory
(next to ``_achilles``) so that auditwheel/delocate can bundle and relocate its
shared-library dependencies. This thin launcher is registered as the ``achilles``
entry point so that ``pip install`` puts the command on the user's PATH.
"""

import os
import subprocess
import sys
from pathlib import Path


def _executable() -> Path:
    exe = Path(__file__).with_name("achilles")
    if os.name == "nt":
        exe = exe.with_suffix(".exe")
    return exe


def main() -> None:
    exe = _executable()
    if not exe.exists():
        sys.exit(
            f"achilles executable not found at {exe}. This wheel may have been "
            "built without the CLI."
        )
    raise SystemExit(subprocess.call([str(exe), *sys.argv[1:]]))
