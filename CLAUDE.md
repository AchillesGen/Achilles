# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the **landing page** for the Achilles neutrino generator project. It is a Sphinx-based static site that serves as a version selector, linking to versioned documentation hosted on GitHub Pages. The landing page lives on the `landing_page` branch (a git worktree of the main Achilles repo).

## Build Commands

```bash
# Install dependencies
pip install -r requirements.txt

# Build the HTML site
sphinx-build -b html . _build/html

# Or equivalently via Make
make html
```

## Architecture

- **conf.py**: Sphinx configuration. Uses `sphinx_book_theme` with extensions `breathe`, `sphinxcontrib.bibtex`, and `sphinx_design`. Loads `versions.json` at build time to populate `html_context` for templates.
- **index.rst**: Single-page landing site with inline HTML/JS that fetches `versions.json` at runtime to populate a version-select dropdown.
- **versions.json**: Generated/maintained externally; lists available documentation versions (e.g., `["latest", "v0.2.0"]`). Used both at build time (conf.py) and at runtime (JS in index.rst).
- **CI**: `.github/workflows/landing.yml` builds on push to `landing_page` and deploys to `gh-pages` branch using `peaceiris/actions-gh-pages` with `keep_files: true` to preserve existing versioned docs.
