#!/usr/bin/env python

from setuptools import setup, find_packages
from os import path
# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='nuChic',
    version='1.0',
    description='Neutrino Event Generator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/jxi24/FNALNeuGen',
    author= 'Joshua Isaacson, \
            William Jay, \
            Alessandro Lovato, \
            Pedro A. Machado, \
            Stefan Prestel, \
            Noemi Rocco, \
            Holger Schulz',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        ],
    packages=find_packages(),
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    install_requires=[
        'numpy',
        'vegas',
        ],
    # Provide executable script to run the main code
    entry_points={
        'console_scripts': [
            'nuChic=nuChic:main',
        ],
    },
)

