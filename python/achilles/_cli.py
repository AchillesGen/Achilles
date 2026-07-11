# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""Console-script launcher for the bundled native ``achilles`` executable.

The compiled C++ ``achilles`` binary is installed inside the package directory
(next to ``_achilles``) so that auditwheel/delocate can bundle and relocate its
shared-library dependencies. This thin launcher is registered as the ``achilles``
entry point so that ``pip install`` puts the command on the user's PATH.

The heavy data set is not shipped in the wheel; it is downloaded on demand into a
user cache (see ``_data.py``). Before launching, we point the binary at that cache via
``ACHILLES_DATA_DIR`` and, if the data is missing, offer to download it first.
"""

import os
import subprocess
import sys
from pathlib import Path

from . import _data


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

    env = dict(os.environ)
    # Resolve data from the shared cache, the same directory the Python module uses.
    env.setdefault("ACHILLES_DATA_DIR", str(_data.cache_dir()))

    # If the data set isn't cached and there's no local ./data, offer to fetch it.
    # download_all() prompts for consent; a decline is non-fatal (the run may still
    # find data via ACHILLES_PATH or the current directory, and if not the binary
    # reports the missing file itself).
    if not _data.is_present() and not (Path.cwd() / "data").is_dir():
        try:
            _data.download_all()
        except RuntimeError as exc:
            print(exc, file=sys.stderr)

    raise SystemExit(subprocess.call([str(exe), *sys.argv[1:]], env=env))
