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
rst_prolog = """
.. include:: <s5defs.rst>
"""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinx_rtd_theme'
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
html_css_files = ['css/s5defs-roles.css']

# -- Breathe configuration ---------------------------------------------------
breathe_default_project = "Achilles"

import subprocess, os

def configureDoxyfile(input_dir, output_dir):
    with open('Doxyfile.in', 'r') as file :
        filedata = file.read()

    filedata = filedata.replace('@DOXYGEN_INPUT_DIR@', input_dir)
    filedata = filedata.replace('@DOXYGEN_OUTPUT_DIR@', output_dir)

    with open('Doxyfile', 'w') as file:
        file.write(filedata)

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

breathe_projects = {}

if read_the_docs_build:
    input_dir = '../include/Achilles'
    output_dir = 'build'
    configureDoxyfile(input_dir, output_dir)
    subprocess.call('doxygen', shell=True)
    breathe_projects['Achilles'] = output_dir + '/xml'

