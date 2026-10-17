.. _FAQ:

##########################
Troubleshooting & FAQ
##########################

.. contents::
   :local:

***************
Build Issues
***************

CMake cannot find HDF5
======================

If CMake reports that HDF5 is not found, ensure that HDF5 is installed and its
location is visible to CMake:

.. code-block:: shell-session

   $ cmake -DHDF5_ROOT=/path/to/hdf5 ..

On systems with multiple HDF5 installations, setting ``HDF5_ROOT`` explicitly
avoids picking up the wrong version.

Missing Fortran compiler
========================

Achilles requires a Fortran compiler for the nuclear model interface. If CMake
reports no Fortran compiler, install ``gfortran`` (or your preferred Fortran
compiler) and ensure it is on your ``PATH``:

.. code-block:: shell-session

   # Ubuntu/Debian
   $ sudo apt install gfortran

   # macOS with Homebrew
   $ brew install gcc

CMake version too old
=====================

Achilles requires CMake 3.12 or newer. Check your version with:

.. code-block:: shell-session

   $ cmake --version

If it is too old, install a newer version from `cmake.org <https://cmake.org/download/>`_
or through your package manager.


*****************
Runtime Errors
*****************

Missing data files
==================

If Achilles exits with an error about missing data files (spectral functions, density
profiles, nucleon configurations, etc.), check that:

1. The paths in your run card are correct and relative to the Achilles data directory.
2. The data files were properly installed. If building from source, data files are
   located in the ``data/`` directory of the source tree and should be accessible from
   the working directory where you run Achilles.

Run card validation errors
==========================

Achilles validates the run card at startup. Common mistakes include:

- **Typos in section names**: YAML keys are case-sensitive. ``NuclearModels`` is not
  the same as ``nuclearmodels``.
- **Incorrect PID values**: Make sure the PDG particle IDs in the ``Processes`` and
  ``Beams`` sections are valid. Common PIDs: electron (11), electron neutrino (12),
  muon neutrino (14).
- **Missing required keys**: Each section has required fields. Refer to the
  :ref:`Run Card Structure` documentation for the complete list.

YAML parsing errors
===================

If you see a YAML parsing error, check for:

- Inconsistent indentation (YAML requires spaces, not tabs).
- Missing colons after keys.
- Unquoted strings that contain special YAML characters (``:  { } [ ] , & * # ? | - < > = ! % @ \``).


*********************
Performance Tips
*********************

Number of events
================

The ``NEvents`` setting controls how many *unweighted* events are produced. The
actual number of phase-space points sampled during integration and unweighting is
much larger. Start with a small number (10--100) to verify your setup, then increase
for production runs.

Unweighting efficiency
======================

If the unweighting efficiency is very low (reported in the log output), events take
longer to generate because many trial points are rejected. Possible remedies:

- **Lower the unweighting percentile**: Setting ``percentile: 95`` instead of ``99``
  discards more high-weight outliers, which raises efficiency at the cost of a small
  bias.
- **Increase the integration accuracy**: A tighter ``Accuracy`` value in the
  ``Initialize`` section produces a better phase-space mapping, improving efficiency.

Parallelism
===========

Achilles does not currently support multi-threaded event generation within a single
run. To parallelise, launch multiple independent runs with different random seeds
and merge the output files afterward.


*****************
Output Issues
*****************

Empty output files
==================

If the output file is empty or contains zero events:

- Check the log for errors during the integration or unweighting phase.
- Ensure that the beam energy is above the kinematic threshold for the requested
  process.
- If hard cuts are enabled, verify that the cut ranges are physically accessible.
  Overly restrictive cuts can reject all events.

Understanding event weights
===========================

In NuHepMC output, the ``"CV"`` weight is the central-value event weight. For
unweighted events, this weight should be close to 1. Weighted events (produced when
unweighting fails for a particular phase-space point) will carry larger weights.

The total cross section can be obtained from the run-level metadata stored in the
HepMC3 ``GenRunInfo`` object. See :ref:`Output Formats and Analysis` for details on
reading the output.


*********************
How-To Guides
*********************

Adding a custom nuclear model
=============================

To implement a custom nuclear model in Fortran:

1. Write a Fortran module that extends the abstract ``model`` type.
2. Register it with the factory in ``nuclear_model_interface.f90``.
3. Add the source file to the CMake build.
4. Reference the model in your run card under ``NuclearModels``.

See the :ref:`Fortran interface <fortran-interface>` documentation for the complete
API reference and a step-by-step example.

Using a custom flux
===================

To use a custom neutrino flux, set the beam type to ``Spectrum`` and provide your
flux histogram file:

.. code-block:: yaml

   Beams:
     - Beam:
         PID: 14
         Beam Params:
           Type: Spectrum
           Histogram: path/to/my_flux.txt

Achilles supports several flux file formats (Achilles, MiniBooNE, T2K) auto-detected
from the file header, as well as HepData YAML files and ROOT histograms (if compiled
with ROOT support). See the :ref:`Run Card Structure` documentation for the full list
of ``Beam Params`` options.
