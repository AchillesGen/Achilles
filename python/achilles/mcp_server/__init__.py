"""Achilles MCP server — Model Context Protocol interface for neutrino event generation."""


def main():
    from .server import main as _main

    _main()


__all__ = ["main"]
