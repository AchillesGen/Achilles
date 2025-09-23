# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Achilles'
copyright = '2025, Achilles Collaboration'
author = 'Joshua Isaacson, William Jay, Alessandro Lovato, Pedro A. Machado, Luke Pickering, Hayden Piwonka, Noemi Rocco, Noah Steinberg'
release = '0.3.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe", "sphinxcontrib.bibtex", "sphinx_design"]

bibtex_bibfiles = ["src/references.bib"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_favicon = '../assets/favicon.ico'
# html_sidebars = {
#     '**': [
#         'about.html',
#         'navigation.html',
#         'relations.html',
#         'searchbox.html',
#         'donate.html',
#     ]
# }
html_theme_options = {
    'logo': {
        "image_light": '../assets/logo.svg',
        "image_dark": '../assets/logo.svg',
    }
}

import json
import os

# Load versions.json
versions_file = os.path.join(os.path.dirname(__file__), "versions.json")
if os.path.exists(versions_file):
    with open(versions_file) as f:
        versions = json.load(f)
else:
    versions = []

# Make it available to templates
html_context = {
    "versions": versions,
}

