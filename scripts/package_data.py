#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""Build the Achilles data-set release assets.

Produces, in the chosen output directory:
  * ``achilles-data.tar.gz``          — the ``data/`` tree, gzip-compressed
  * ``achilles-data.tar.gz.sha256``   — its SHA-256 checksum (``sha256sum`` format)

Upload both as assets on the GitHub Release whose tag matches ``DATA_RELEASE`` in
``CMakeLists.txt`` (e.g. ``achilles-data-v1``). The runtime downloader
(``python/achilles/_data.py``) fetches them from that release and verifies the
checksum before extracting.

Usage:
    python scripts/package_data.py [--data-dir data] [--out-dir dist]
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import tarfile
from pathlib import Path

ASSET_NAME = "achilles-data.tar.gz"


def build(data_dir: Path, out_dir: Path) -> Path:
    if not data_dir.is_dir():
        sys.exit(f"data directory not found: {data_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    tarball = out_dir / ASSET_NAME

    # arcname="data" so the archive extracts to <cache>/data/..., matching how run
    # cards reference files (e.g. "data/configurations/...").
    print(f"Packing {data_dir} -> {tarball}")
    with tarfile.open(tarball, "w:gz") as tar:
        tar.add(data_dir, arcname="data")

    digest = hashlib.sha256(tarball.read_bytes()).hexdigest()
    sha_file = out_dir / f"{ASSET_NAME}.sha256"
    sha_file.write_text(f"{digest}  {ASSET_NAME}\n")
    print(f"SHA-256: {digest}")
    print(f"Wrote {tarball} ({tarball.stat().st_size / 1e6:.1f} MB) and {sha_file}")
    return tarball


def main(argv: list[str] | None = None) -> int:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description="Build Achilles data release assets.")
    parser.add_argument("--data-dir", type=Path, default=repo_root / "data")
    parser.add_argument("--out-dir", type=Path, default=repo_root / "dist")
    args = parser.parse_args(argv)
    build(args.data_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
