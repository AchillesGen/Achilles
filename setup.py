#!/usr/bin/env python
""" Setup script for the nuchic neutrino event generator."""

import io
import os
from os.path import join
from os.path import dirname
import re
import sys
import platform
import subprocess

from distutils.version import LooseVersion
import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


def read(*names, **kwargs):
    """ Read files. """
    with io.open(
            join(dirname(__file__), *names),
            encoding=kwargs.get('encoding', 'utf8')
    ) as fhs:
        return fhs.read()


class CMakeExtension(Extension):
    """ Create a CMake extension type for python. """
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """ Build a CMake project during the building of the python library. """
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == 'Windows':
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name='nuchic',
    version='0.0.1',
    description='Neutrino Event Generator',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    url='https://github.com/jxi24/FNALNeuGen',
    author='Joshua Isaacson, \
            William Jay, \
            Alessandro Lovato, \
            Pedro A. Machado, \
            Noemi Rocco,',
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
    include_package_data=True,
    python_requires='>=3.5, <4',
    install_requires=[
        'numpy',
        'vegas',
        'pandas',
        'matplotlib',
        'tqdm',
        'pyyaml',
        'scipy',
        'h5py',
    ],
    entry_points={'console_scripts': ['nuchic = nuchic.main:nuchic', ], },
    ext_modules=[CMakeExtension('nuchic')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    package_data={'': ['data/*', 'data/qe/*', 'pke/*', 'configurations/*', 'template.yml']},
    extras_require={
        'test': ['pytest', 'coverage', 'pytest-cov'],
    },
)
