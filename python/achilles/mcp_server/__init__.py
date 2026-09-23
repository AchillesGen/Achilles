# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""Achilles MCP server — Model Context Protocol interface for neutrino event generation."""


def main():
    from .server import main as _main

    _main()


__all__ = ["main"]
