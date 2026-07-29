# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""On-demand download of the Achilles data set.

The wheel ships no heavy data (the ``data/`` tree is 227 MB, past PyPI's 100 MB
per-file limit). Instead the full data set is published as a GitHub Release asset and
downloaded, **only with explicit user consent**, into a cross-platform user cache
directory. Both the Python module and the bundled ``achilles`` CLI find it there via
the ``ACHILLES_DATA_DIR`` environment variable, which is added to the C++ search path
in ``Filesystem::AchillesPath`` (``src/Achilles/System.cc``).

The cache directory is chosen by :func:`platformdirs.user_data_dir`, so there is no
OS-dependent path logic in this file.

Data releases are versioned **independently of the code**. The required achilles-data
tag is set once in ``CMakeLists.txt`` (``ACHILLES_DATA_RELEASE``) and baked into the
generated ``_data_release.py`` at configure time. Data changes rarely, so most code
releases reuse the same pin.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

import platformdirs

try:
    # Generated at CMake configure time; canonical value lives in CMakeLists.txt.
    from ._data_release import DATA_RELEASE as _DEFAULT_DATA_RELEASE
except ImportError:
    # Raw source tree that was never configured by CMake. Keep this fallback in sync
    # with the ACHILLES_DATA_RELEASE default in CMakeLists.txt.
    _DEFAULT_DATA_RELEASE = "achilles-data-v1"

_REPO = "AchillesGen/Achilles"

_ENV_DATA_DIR = "ACHILLES_DATA_DIR"
_ENV_DATA_RELEASE = "ACHILLES_DATA_RELEASE"
#: Full override for the tarball URL (used by tests / mirrors).
_ENV_DATA_URL = "ACHILLES_DATA_URL"

_TARBALL_NAME = "achilles-data.tar.gz"


def data_release() -> str:
    """The achilles-data tag to fetch (env-overridable, else the configured pin)."""
    return os.environ.get(_ENV_DATA_RELEASE, _DEFAULT_DATA_RELEASE)


def cache_dir() -> Path:
    """Root of the data cache. Contains the ``data/`` subtree once populated.

    Honors ``ACHILLES_DATA_DIR`` if the user set it, otherwise a platform-appropriate
    user data directory. Never touches the filesystem.
    """
    override = os.environ.get(_ENV_DATA_DIR)
    if override:
        return Path(override)
    return Path(platformdirs.user_data_dir("achilles", "AchillesGen"))


def data_root() -> Path:
    """The ``data/`` directory inside the cache (what files are referenced against)."""
    return cache_dir() / "data"


def _stamp_path() -> Path:
    # Keyed on the data release so bumping the pin forces a fresh download.
    return cache_dir() / f".achilles-{data_release()}"


def is_present() -> bool:
    """True if the pinned data release has been fully downloaded."""
    return _stamp_path().is_file() and data_root().is_dir()


def _tarball_url() -> str:
    override = os.environ.get(_ENV_DATA_URL)
    if override:
        return override
    return (
        f"https://github.com/{_REPO}/releases/download/"
        f"{data_release()}/{_TARBALL_NAME}"
    )


def _sha256_url() -> str:
    return _tarball_url() + ".sha256"


def _fetch_expected_sha256() -> str | None:
    """Download the ``.sha256`` sidecar for the tarball, or None if unavailable."""
    try:
        with urllib.request.urlopen(_sha256_url()) as resp:  # noqa: S310 (trusted host)
            first = resp.read().decode().split()
            return first[0] if first else None
    except (urllib.error.URLError, OSError):
        return None


def _remote_size() -> int | None:
    try:
        req = urllib.request.Request(_tarball_url(), method="HEAD")
        with urllib.request.urlopen(req) as resp:  # noqa: S310
            length = resp.headers.get("Content-Length")
            return int(length) if length else None
    except (urllib.error.URLError, OSError, ValueError):
        return None


def _human_size(num: int | None) -> str:
    if num is None:
        return "an unknown amount of data"
    size = float(num)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}"
        size /= 1024
    return f"{size:.0f} GB"


def _confirm(prompt: str) -> bool:
    if not sys.stdin or not sys.stdin.isatty():
        # Non-interactive: never assume yes. Caller must pass assume_yes explicitly.
        return False
    try:
        return input(f"{prompt} [y/N] ").strip().lower() in ("y", "yes")
    except EOFError:
        return False


def _safe_extract(tar: tarfile.TarFile, dest: Path) -> None:
    """Extract guarding against path-traversal (``..``/absolute) members.

    Uses the stdlib ``data`` filter on 3.12+, and an equivalent manual check on the
    3.10/3.11 minimum we support.
    """
    dest = dest.resolve()
    if hasattr(tarfile, "data_filter"):
        tar.extractall(dest, filter="data")
        return
    for member in tar.getmembers():
        target = (dest / member.name).resolve()
        if not str(target).startswith(str(dest)):
            raise RuntimeError(f"Refusing to extract unsafe path: {member.name!r}")
    tar.extractall(dest)  # noqa: S202 (members validated above)


def download_all(*, assume_yes: bool = False) -> Path:
    """Download and unpack the pinned data release into the cache. Returns the cache dir.

    Raises RuntimeError on verification failure or if the user declines.
    """
    if is_present():
        print(f"Achilles data ({data_release()}) already present in {cache_dir()}")
        return cache_dir()

    size = _remote_size()
    url = _tarball_url()
    if not assume_yes:
        prompt = (
            f"Download the Achilles data set '{data_release()}' "
            f"({_human_size(size)}) from\n  {url}\ninto {cache_dir()}?"
        )
        if not _confirm(prompt):
            raise RuntimeError(
                "Download declined. Run `achilles-data get` to fetch the data set."
            )

    expected = _fetch_expected_sha256()
    dest = cache_dir()
    dest.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(dir=dest, suffix=".tar.gz.part", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        try:
            print(f"Downloading {url} ...")
            digest = hashlib.sha256()
            with urllib.request.urlopen(url) as resp:  # noqa: S310
                while chunk := resp.read(1 << 20):
                    tmp.write(chunk)
                    digest.update(chunk)
            tmp.flush()

            actual = digest.hexdigest()
            if expected is None:
                print(
                    "Warning: no published checksum found; skipping integrity check.",
                    file=sys.stderr,
                )
            elif actual != expected:
                raise RuntimeError(
                    f"Checksum mismatch for {url}\n  expected {expected}\n  got      {actual}"
                )

            print(f"Extracting into {dest} ...")
            with tarfile.open(tmp_path, "r:gz") as tar:
                _safe_extract(tar, dest)
        finally:
            tmp_path.unlink(missing_ok=True)

    _stamp_path().write_text(expected or "unverified")
    print(f"Achilles data ready in {dest}")
    return dest


def status() -> None:
    present = is_present()
    print(f"Data release : {data_release()}")
    print(f"Cache dir    : {cache_dir()}")
    print(f"Status       : {'present' if present else 'MISSING'}")
    if not present:
        print("Run `achilles-data get` to download the data set.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="achilles-data",
        description="Manage the on-demand Achilles data set.",
    )
    parser.add_argument(
        "-y", "--yes", action="store_true", help="skip the confirmation prompt"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status", help="show whether the data set is downloaded")
    sub.add_parser("get", help="download the pinned data set")
    sub.add_parser("path", help="print the data cache directory")

    args = parser.parse_args(argv)
    if args.command == "status":
        status()
        return 0
    if args.command == "path":
        print(cache_dir())
        return 0
    if args.command == "get":
        try:
            download_all(assume_yes=args.yes)
        except RuntimeError as exc:
            print(exc, file=sys.stderr)
            return 1
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
