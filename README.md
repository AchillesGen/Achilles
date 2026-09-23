<h1 align="center">
  <img src="/assets/logo.svg" alt="Achilles" width="25%"/>
</h1>

[![CMake Build Matrix](https://github.com/AchillesGen/Achilles/actions/workflows/cmake.yml/badge.svg)](https://github.com/AchillesGen/Achilles/actions/workflows/cmake.yml)

[![codecov](https://codecov.io/gh/Achilles/AchillesGen/branch/main/graph/badge.svg?token=Xq2sJ4kv5L)](https://codecov.io/gh/AchillesGen/Achilles)

## Introduction

Achilles (**A** **CHI**cago **L**and **L**epton **E**vent **S**imulator) is a Monte Carlo event generator for the simulation of lepton-nucleus and neutrino-nucleus interactions.

The main goals of the Achilles event generator are:

1. Theory driven: ensure that the guiding principles are to reproduce our understanding of the Standard Model with appropriate uncertainties
2. Leverage experiences from the LHC event generator community
3. Develop a modular and extensible neutrino event generator
4. Provide automated Beyond the Standard Model calculations for neutrino experiments
5. Fully evaluate theory uncertainties

For comprehensive documentation, please visit the [Achilles documentation site](https://achillesgen.github.io/Achilles/).
The full manual can be found [here](https://achillesgen.github.io/Achilles/latest/).

## Building Achilles

The Achilles code uses CMake to provide a platform-agnostic installation procedure. A C++17 compiler, a Fortran90 compiler (we recommend gfortran), and [HDF5](https://www.hdfgroup.org/solutions/hdf5/) are required. All additional dependencies are automatically downloaded and built by CMake.

To build Achilles with the default options:

```bash
cmake -S . -B build
cmake --build build -- -jN
```

To install Achilles locally (no superuser permissions required):

```bash
cmake --install build --prefix build/install
```

The install location can be changed by passing a different path to `--prefix` or by setting `CMAKE_INSTALL_PREFIX` at configure time.

For the full list of build options and optional dependencies (including the Sherpa interface for BSM processes and tau decays), see the [Getting Started](https://achillesgen.github.io/Achilles/latest/src/getting-started.html) documentation.

### CMake Options

A summary of commonly used CMake options:

| Option                           | Meaning                                                                          | Default |
| -------------------------------- | -------------------------------------------------------------------------------- | ------- |
| `ACHILLES_ENABLE_TESTING`        | Build the Achilles test suite                                                    | OFF     |
| `ACHILLES_ENABLE_CASCADE_TEST`   | Build the executable to only run the cascade (hadron-nucleus cross section or transparency) | OFF     |
| `ACHILLES_ENABLE_POTENTIAL_TEST` | Build executable to test different potentials                                    | OFF     |
| `ACHILLES_ENABLE_SHERPA`         | Build the Sherpa interface for BSM and tau decays                                | OFF     |
| `ACHILLES_BUILD_DOCS`            | Build the Achilles manual                                                        | OFF     |
| `CMAKE_BUILD_TYPE`               | Build type: Release, Debug, or RelWithDebInfo                                    | Release |

## Running Achilles

The main Achilles executable can be found at `bin/achilles` after building. Running with no arguments looks for a `run.yml` run card in the current directory. To use a different run card:

```bash
./bin/achilles <run_card>
```

For details on the run card structure, runtime options, output formats, and Docker support, see the [full documentation](https://achillesgen.github.io/Achilles/latest/src/getting-started.html).

## Testing Achilles

If the code was configured with `ACHILLES_ENABLE_TESTING=ON`, the test suite can be run using:

```
./test/achilles-testsuite
```

## Citing Achilles

If you use Achilles, please cite:

```
@article{Isaacson:2022cwh,
    author = "Isaacson, Joshua and Jay, William I. and Lovato, Alessandro and Machado, Pedro A. N. and Rocco, Noemi",
    title = "{ACHILLES: A novel event generator for electron- and neutrino-nucleus scattering}",
    eprint = "2205.06378",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FERMILAB-PUB-22-411-T, MIT-CTP/5428",
    month = "5",
    year = "2022"
}
```

If you use Achilles for a BSM calculation, please additionally cite:

```
@article{Isaacson:2021xty,
    author = {Isaacson, Joshua and H\"oche, Stefan and Lopez Gutierrez, Diego and Rocco, Noemi},
    title = "{Novel event generator for the automated simulation of neutrino scattering}",
    eprint = "2110.15319",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FERMILAB-PUB-21-537-T, MCNET-21-31",
    doi = "10.1103/PhysRevD.105.096006",
    journal = "Phys. Rev. D",
    volume = "105",
    number = "9",
    pages = "096006",
    year = "2022"
}
```

```
@article{Hoche:2014kca,
    author = {H\"oche, Stefan and Kuttimalai, Silvan and Schumann, Steffen and Siegert, Frank},
    title = "{Beyond Standard Model calculations with Sherpa}",
    eprint = "1412.6478",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "SLAC-PUB-16170, IPPP-14-105, DCPT-14-210, MCNET-14-35",
    doi = "10.1140/epjc/s10052-015-3338-4",
    journal = "Eur. Phys. J. C",
    volume = "75",
    number = "3",
    pages = "135",
    year = "2015"
}
```

```
@article{Gleisberg:2008fv,
    author = "Gleisberg, Tanju and Hoeche, Stefan",
    title = "{Comix, a new matrix element generator}",
    eprint = "0808.3674",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "SLAC-PUB-13232, IPPP-08-31, DCPT-08-62, MCNET-08-08",
    doi = "10.1088/1126-6708/2008/12/039",
    journal = "JHEP",
    volume = "12",
    pages = "039",
    year = "2008"
}
```
