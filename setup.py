#!/usr/bin/env python
""" Setup script for the nuchic neutrino event generator. """

# from setuptools import setup, find_packages, Extension
import io
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext
import setuptools
from numpy.distutils.core import setup, Extension


def read(*names, **kwargs):
    """ Read files. """
    with io.open(
            join(dirname(__file__), *names),
            encoding=kwargs.get('encoding', 'utf8')
    ) as fhs:
        return fhs.read()


EXT1 = Extension(name='xsec', sources=['src/nuchic/xsec.pyf',
                                       'src/nuchic/currents_opt_v1.f90',
                                       'src/nuchic/mathtool.f90',
                                       'src/nuchic/nform.f90',
                                       'src/nuchic/xsec.f90'])

setup(
    name='nuchic',
    version='1.0',
    description='Neutrino Event Generator',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    url='https://github.com/jxi24/FNALNeuGen',
    author='Joshua Isaacson, \
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
        'Programming Language :: Python :: 3',
    ],
    packages=setuptools.find_packages(where='src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    python_requires='>=3.5, <4',
    install_requires=[
        'numpy',
        'vegas',
        'h5py',
        'scipy',
        'pandas',
        'sklearn',
        'absl-py',
        'matplotlib',
        'tqdm',
        'pyyaml',
    ],
    # Provide executable script to run the main code
    entry_points={'console_scripts': ['nuchic = nuchic.main:nu_chic', ], },
    ext_modules=[EXT1],
    package_data={'': ['data/*', 'data/qe/*', 'pke/*', 'configurations/*', 'template.yml']},
    extras_require={
        'test': ['pytest', 'coverage', 'pytest-cov'],
    },
)
