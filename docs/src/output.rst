.. _Output Formats and Analysis:

############################
Output Formats & Analysis
############################

.. contents::
   :local:

Achilles writes generated events to disk in one of three formats, configured via the
``Format`` key in the ``Output`` sub-node of the run card (see :ref:`Run Card Structure`
for details). This page describes each format and how to use the output with external
analysis tools.

*****************************
NuHepMC Format (Recommended)
*****************************

The recommended output format is `NuHepMC <https://github.com/NuHepMC/Spec>`_, a set of
conventions built on top of the `HepMC3 <https://gitlab.cern.ch/hepmc/HepMC3>`_ event
record. NuHepMC files are readable by any HepMC3-compatible tool while also carrying
neutrino-experiment-specific metadata.

Achilles sets the following NuHepMC signal conventions:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Convention
     - Description
   * - ``G.C.4``
     - Cross-section units and target scaling are specified in the run info.
   * - ``G.C.6``
     - Particle and vertex status codes are defined in the run info.
   * - ``E.C.1``
     - Each event carries a process ID identifying the interaction channel.
   * - ``E.C.4``
     - Each event records the incoming neutrino flux value.
   * - ``E.C.5``
     - The lab-frame interaction position is stored on each event.
   * - ``V.C.1``
     - The primary interaction vertex is identified by its status code.

Cross sections are stored in picobarns per atom. Event weights are accessible through the
``"CV"`` weight name. The output file can optionally be gzip-compressed by setting
``Zipped: true`` in the run card.

Reading NuHepMC files
=====================

Any HepMC3 reader can open NuHepMC output. In Python:

.. code-block:: python

   import pyHepMC3 as hep

   reader = hep.ReaderAscii("achilles.hepmc")
   evt = hep.GenEvent()
   while not reader.failed():
       reader.read_event(evt)
       if reader.failed():
           break
       # Access particles, vertices, weights, etc.
       print(evt.weight("CV"))

The `NuHepMC C++ library <https://github.com/NuHepMC/cpputils>`_ provides additional
helper functions for extracting convention-specific metadata.

****************************
Using Output with NUISANCE
****************************

`NUISANCE <https://github.com/NUISANCEMC/nuisance>`_ is a cross-section comparison
framework for neutrino event generators. It can read Achilles NuHepMC output directly
and compare generated distributions against published experimental data.

`NUISANCE3 <https://github.com/NUISANCEMC/nuisance3>`_ is the next-generation rewrite
of NUISANCE with native NuHepMC support and a modernised analysis framework. Both
versions accept HepMC3/NuHepMC files as input.

A typical workflow:

1. Generate events with Achilles in NuHepMC format.
2. Pass the output file to NUISANCE to produce comparison plots against data.
3. Iterate on the run card or nuclear model settings as needed.

***********************
HepMC Format (Legacy)
***********************

.. deprecated:: 0.2
   The plain HepMC format is retained for backward compatibility. New users should use
   NuHepMC instead.

Plain HepMC3 output without NuHepMC metadata. The event record structure is identical but
the neutrino-specific convention attributes are omitted. Select this format by setting
``Format: HepMC`` in the run card.

****************************
Achilles Native Format
****************************

The native Achilles binary format is a compact representation intended for internal
debugging and post-processing within the Achilles framework. It is not designed for
external analysis and its layout may change between versions. Select this format by
setting ``Format: Achilles`` in the run card.
