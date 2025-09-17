<h1 align="center">
  <img src="/assets/logo.svg" alt="Achilles" width="25%"/>
</h1>

[![CMake Build Matrix](https://github.com/AchillesGen/Achilles/actions/workflows/cmake.yml/badge.svg)](https://github.com/AchillesGen/Achilles/actions/workflows/cmake.yml)

[![codecov](https://codecov.io/gh/Achilles/AchillesGen/branch/main/graph/badge.svg?token=Xq2sJ4kv5L)](https://codecov.io/gh/AchillesGen/Achilles)

## Introduction

Achilles (A CHIcago Land Lepton Event Simulator) is a modern theory driven lepton event generator.
The focus of the generator is to simulate electron-nucleus and neutrino-nucleus scattering.
The design of the code is based on the following principles:
1. Modular framework to switch in different models
2. Easy extension by the users
3. Theory driven with appropriate uncertainties
4. Provide automated BSM calculations for neutrino experiments

Additional details can be found in the Achilles [wiki](https://github.com/AchillesGen/Achilles/wiki).

## Why a new generator?

TODO: Add details in this section

## Building Achilles

In this section the basic method of building the Achilles code is provided.
For further details and options, please refer to [build details](https://github.com/AchillesGen/Achilles/wiki/Build-Details).
The Achilles code uses CMake as a means to provide a platform agnositic installation procedure.

The default options for the building of Achilles requires HepMC3
and Sherpa. The HepMC3 code provides a means to output events in the convention dictated by the [NuHepMC](https://github.com/NuHepMC/Spec) standard.
The Sherpa interface allows for the simulation of beyond the Standard Model (BSM) processes. Details on obtaining
these codes can be found in the next [section](#-optional-dependencies).

To build Achilles with these default options can be done with:
```bash
mkdir build && cd build
cmake .. -DACHILLES_ENABLE_SHERPA=ON -DSHERPA_ROOT_DIR=/path/to/sherpa
make -jN
```

If the HepMC3 cmake files are not within the CMake module path, you can add the `-DHepMC3_DIR=/path/to/hepmc3/cmake/files`
to the above `cmake` command. If HepMC3 is not present, Achilles will install HepMC3 for you. Additional details and optional dependencies can be found below.

### Optional Dependencies

#### HepMC3

By default, Achilles will install HepMC3 for you.
Alternatively, you can install it yourself.

The HepMC3 code can be found [here](https://gitlab.cern.ch/hepmc/HepMC3), and has details on building and
installing the code. Achilles requires HepMC3 version 3.2.5 or newer.

HepMC3 provides a C++ and python interface for writing HepMC3 files based on the [arxiv:1912.08005](https://arxiv.org/abs/1912.08005).
The HepMC3 is supported and maintained by the LHC and heavy ion communities. This has become a
standard in the HEP event generator community.

For details on the additions to the HepMC3 standard for colliders to neutrino physics see [here](https://github.com/NuHepMC/Spec).

To disable the requirement of HepMC3, add the option `-DENABLE_HEPMC3=OFF` to the cmake command.

#### Sherpa

The leptonic currents are calculated as described in [arxiv:2110.15319](https://arxiv.org/abs/2110.15319). This involves calculating
temrs using the Berends-Giele recursion relations in arbitrary models. The calculation of these is
implemented into the Comix matrix element generator within the Sherpa codebase.

The required version of Sherpa is in the process of being made public, but can be supplied upon request to the
Achilles authors.
Note that to enable UFO support from Sherpa, add the option `--enable-ufo' to the configure command.

To disable the requirement of Sherpa, add the option `-DENABLE_BSM=OFF` to the cmake command.

### CMake Options

A non-exhaustive list of useful CMake options includes:

| Option                  | Meaning                                                                         |
| ------                  | -------                                                                         |
| `ACHILLES_ENABLE_TESTING`        | Build the Achilles test suite                                                   |
| `ACHILLES_ENABLE_GZIP`           | Compile the code with the ability to directly compress event files              |
| `ACHILLES_ENABLE_CASCADE_TEST`   | Build the executable to only run the cascade (pA cross section or transparency) |
| `ACHILLES_ENABLE_POTENTIAL_TEST` | Build executable to test different potentials                                   |
| `CMAKE_BUILD_TYPE=Debug`         | Build the executable with debug symbols enabled                                |
## Running Achilles

The main Achilles executable can be found at `bin/achilles` after building the code. Running `./bin/achilles --help` will provide all the different command line options available to the user. Currently, these are:

```
    Usage:
      achilles [<input>] [-v | -vv] [-s | --sherpa=<sherpa>...]
      achilles --display-cuts
      achilles --display-ps
      achilles --display-ff
      achilles --display-int-models
      achilles --display-nuc-models
      achilles (-h | --help)
      achilles --version

    Options:
      -v[v]                                 Increase verbosity level.
      -h --help                             Show this screen.
      --version                             Show version.
      -s <sherpa> --sherpa=<sherpa>         Define Sherpa option.
      --display-cuts                        Display the available cuts
      --display-ps                          Display the available phase spaces
      --display-ff                          Display the available form factors
      --display-int-models                  Display the available cascade interaction models
      --display-nuc-models                  Display the available nuclear interaction models
```

The options `--display-cuts`, `--display-ps`, and `--display-ff` will output the available options for
each case and then exit the code. For example, running `./bin/achilles --display-cuts` produces the
following output (splash screen suppressed for brevity):

```
Registered Single Particle cuts:
  - AngleTheta
  - ETheta2
  - Energy
  - Momentum
  - TransverseMomentum
Registered Two Particle cuts:
  - DeltaTheta
  - InvariantMass
```

These options for different cuts can be expressed in the run card as described [below](#-run-card), and
in more details in the [wiki](https://github.com/AchillesGen/Achilles/wiki) and the manual.

### Runtime Options

#### Run card

The run card consists of nine major sections describing how the generation is to be carried out.
These sections are:
1. The main event section
2. The process section
3. The initialization of the random number generator and precision of the integrator section
4. The unweighting method to use
5. The incoming beam
6. Settings for the cascade
7. Settings for the nuclear interaction model
8. Settings for the nucleus
9. Any cuts to apply during the generation of the events

Each of these sections are described below and in greater detail in the
[wiki](https://github.com/AchillesGen/Achilles/wiki).

The _Main_ section contains options:
 - The number of events (`NEvents`)
 - If cuts should be applied at the generation level (`HardCuts`)
 - The output (`Output`), which contains sub-options:
    - The event output format (`Format`, currently options are "HepMC3" and "Achilles")
    - The name of the output file (`Name`)
    - If the file should be written as a gzip file or not (`Zipped`)

The _Process_ section contains information needed to generate the leptonic current for a given physics model.
This contains the options for:
 - The physics model (`Model`)
 - The output leptonic states as a list of particle IDs (`Final States`)

The _Initialization_ section describes the initialization of the generator, and contains:
 - The random seed to use for event generation for reproducibility (`Seed`)
 - The accuracy for the warm-up run of the integrator to achieve before generating events (`Accuracy`)

The _Unweighting_ section sets up the methodology for unweighting the events. This has one required setting
as the `Name` of the unweighting procedure. Each unweighting procedure has their own set of options
described in detail in the [wiki](https://github.com/AchillesGen/Achilles/wiki/Unweighting).

The _Beams_ section provides the means to setup all possible incoming neutrino fluxes.
Currently, only a single flavor incoming beam is supported. The options available for the beam
depends on the type of beam and are explained in detail
in the [wiki](https://github.com/AchillesGen/Achilles/wiki/Beams).

The _Cascade_ section determines the setup of the cascade. The options used to define the cascade are:
 - If the cascade should be ran (`Run`)
 - A sub-section on the calculation of particle interactions to use. This requires the `Name` of the
   interaction model, which can be found using `./bin/achilles --display-int-models`. Additional details
   for the settings for each model can be found in
   the [wiki](https://github.com/AchillesGen/Achilles/wiki/Cascade).
 - The maximum step size to take during the cascade
 - The probability model for determining interactions.
   Currently, only `Cylinder` and `Gaussian` are implemented.
 - If the nucleons should be propagated in a nuclear potential (`PotentialProp`)

The next section is the _Nuclear Model_ section. Here the definition of the nuclear model used for the
primary interaction is defined. The required options are:
 - The model name (`Model`)
 - The file to load the form factors from (`FormFactorFile`). Details of this file can be found in the following
   section.
 - Additional required options depend on the nuclear model used
   and can be found in the [wiki](https://github.com/AchillesGen/Achilles/wiki/Nuclear-Models).

The _Nucleus_ section defines the nucleus for interactions. Currently, only a single isotope and nucleus is
supported to be run at a time. The required options are:
 - The name of the nucleus given as the number of nucleons followed by the chemical symbol (_i.e._ "12C").
 - The Fermi momentum is needed.
 - The setup for the density and configuration.
   Details can be found in the [wiki](https://github.com/AchillesGen/Achilles/wiki/Nucleus).
 - The Fermi gas mode for the cascade. Current options are "Local" and "Global".
 - The nuclear potential to use.
   Details can be found in the [wiki](https://github.com/AchillesGen/Achilles/wiki/Nucleus).

The last section is the _Hard Cuts_ section and defines the cuts to be made on the particles after the
generation of the phase space, but before the cascade. These are used for example to limit the phase
space generated for electron scattering experiments like e4v to more efficiently generate events.
The details of this section are laid out in the [wiki](https://github.com/AchillesGen/Achilles/wiki/Hard-Cuts).

#### Form factors

The form factor file contains the list of the form factors to use, and the parameters for the different
parameterization. Currently, the form factors implemented are:
 - Vector:
    - Dipole
    - Kelly
    - BBBA
    - ArringtonHill
 - Axial:
    - Dipole
 - Coherent:
    - Helm
    - Lovato (Carbon only)

For additional details on the parameters for each form factor, see the [wiki](https://github.com/AchillesGen/Achilles/wiki/Form-Factors).

### Adding models to Achilles (via Sherpa)

The Beyond the Standard Model handling within Achilles is handled via an interface to Sherpa and Comix.
Therefore, in order to add a model to Achilles, you have to process the UFO files through the Sherpa interface.
This can be done with the command `Sherpa-generate-model`, which takes as an input the path to a UFO model
file. Additionally, the model needs to include modifications to handle the interactions with the nucleus which
are currently not automated by FeynRules. Further details can be found in the [wiki](https://github.com/AchillesGen/Achilles/wiki/BSM).

The UFO files for the Dark Neutrino portal model () are included in the repository in the folder `UFO`.
To add this model to be available to Achilles, run the command `Sherpa-generate-model --ncore=N UFO/DarkNeutrinoPortal_Dirac_UFO`. An example run card and parameter card are also provided as `run_hnl.yml` and `hnl_parameters.dat`. Events can be generated with this example file using `./bin/achilles run_hnl.yml`.

## Testing Achilles

If the code was configured with `ACHILLES_ENABLE_TESTING=ON`, the test suite can be run using

```
./test/achilles-testsuite
```
after compilation.
Tests for the cascade and for potential potential are toggled on or off at the time of compilation (see [CMake Options](#cmake-options) above.)

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

If you use Achilles for a BSM calculation, please cite the following three references:

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
