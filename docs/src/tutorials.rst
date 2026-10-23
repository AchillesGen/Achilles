.. _Tutorials:

#########
Tutorials
#########

.. contents::
   :local:

These tutorials walk through running Achilles with the example run cards shipped in
the ``examples/`` directory. Each tutorial explains the run card configuration, shows
how to launch the simulation, and describes the expected output.

Before starting, make sure Achilles is built and the executable is available (see
:ref:`Building Achilles`).


*********************************************************
Tutorial 1: Quasielastic Neutrino Scattering on Carbon
*********************************************************

This tutorial uses the run card ``examples/run_nue_12C_SF_1000.yml``, which simulates
charged-current quasielastic (CCQE) electron-neutrino scattering on a :sup:`12`\ C
target at 1 GeV using a spectral function nuclear model.

The run card
============

.. code-block:: yaml

   Main:
     NEvents: 10
     HardCuts: false
     RunDecays: false
     Output:
         Format: NuHepMC
         Name: achilles_nue_12C_SF_1000.hepmc
         Zipped: True

- **NEvents**: Number of unweighted events to generate. Increase this for production
  runs (e.g. 50000).
- **HardCuts**: No kinematic cuts are applied at the hard scattering level.
- **Output**: Events are written in NuHepMC format, gzip-compressed.

.. code-block:: yaml

   Processes:
     - Leptons: [12, [11]]

The process line defines the leptonic current: an incoming :math:`\nu_e` (PID 12)
producing an outgoing :math:`e^-` (PID 11). This is charged-current scattering.

.. code-block:: yaml

   Beams:
     - Beam:
         PID: 12
         Beam Params:
           Type: Monochromatic
           Energy: 1000

A monochromatic :math:`\nu_e` beam at 1000 MeV.

.. code-block:: yaml

   Cascade:
     Run: False

The intranuclear cascade is disabled. The output contains only the primary hard
scattering products. See :ref:`Tutorial 3 <tutorial-cascade>` for enabling the cascade.

.. code-block:: yaml

   NuclearModels:
     - NuclearModel:
         Model: QESpectral
         FormFactorFile: "FormFactors.yml"
         SpectralP: data/Spectral_Functions/pke12p_tot.data
         SpectralN: data/Spectral_Functions/pke12n_tot.data
         Ward: None

The nuclear model is ``QESpectral`` -- quasielastic scattering computed with proton
and neutron spectral functions for :sup:`12`\ C. No Ward gauge correction is applied.

.. code-block:: yaml

   Nuclei:
     - Nucleus:
         Name: 12C
         Binding: 8.6
         Fermi Momentum: 225
         ...

The target nucleus is :sup:`12`\ C with a binding energy of 8.6 MeV and Fermi
momentum of 225 MeV. The density profile and nucleon configurations are loaded from
the data files.

Running the simulation
======================

From the build directory:

.. code-block:: shell-session

   $ ./achilles ../examples/run_nue_12C_SF_1000.yml

Achilles will:

1. Initialise the adaptive Monte Carlo integrator.
2. Run an optimisation phase to map out the phase space.
3. Generate the requested number of unweighted events.
4. Write the output to ``achilles_nue_12C_SF_1000.hepmc.gz``.

Progress is logged to the terminal and to ``achilles.log``.

Expected output
===============

The output file ``achilles_nue_12C_SF_1000.hepmc.gz`` contains unweighted events
in NuHepMC format. Each event has:

- An incoming :math:`\nu_e` and a target nucleus (:sup:`12`\ C).
- A primary vertex with an outgoing :math:`e^-` and a knocked-out nucleon.
- Spectator nucleons from the nuclear configuration.
- Event weights and cross-section metadata.


*************************************************
Tutorial 2: Electron Scattering on Hydrogen
*************************************************

This tutorial uses ``examples/run_electron_1H_SF_1000.yml``, which simulates
electron scattering on a free proton (hydrogen) at 1 GeV. This is a simpler
configuration that serves as a good starting point for understanding the code.

Key differences from Tutorial 1
================================

.. code-block:: yaml

   Processes:
     - Leptons: [11, [11]]

Both the incoming and outgoing leptons are electrons (PID 11). This is
electromagnetic (electron) scattering rather than neutrino scattering.

.. code-block:: yaml

   Main:
     HardCuts: true

   HardCuts:
     - Type: AngleTheta
       PIDs: 11
       range: [10, 90]

Hard cuts are enabled. The outgoing electron is required to have a polar angle
:math:`\theta` between 10 and 90 degrees. This mimics a typical electron scattering
experiment acceptance.

.. code-block:: yaml

   Nuclei:
     - Nucleus:
         Name: 1H

The target is hydrogen (a free proton). Because there is no nuclear medium, the
spectral function and Fermi gas settings have no effect -- Achilles automatically
treats hydrogen as a single-nucleon target.

Running the simulation
======================

.. code-block:: shell-session

   $ ./achilles ../examples/run_electron_1H_SF_1000.yml

Expected output
===============

The output file ``achilles_electron_1000_SF_1H.hepmc.gz`` contains events with:

- An incoming and outgoing electron.
- A recoil proton.
- No nuclear effects (no spectator nucleons, no cascade).

This provides a clean baseline for validating form factor implementations against
known electron-proton scattering data.


.. _tutorial-cascade:

*************************************
Tutorial 3: Enabling the Cascade
*************************************

The intranuclear cascade (INC) propagates hadronic final-state particles through the
nuclear medium, allowing them to rescatter, produce pions, or be absorbed. This
tutorial shows how to modify a run card to enable the cascade.

Modifying the run card
======================

Starting from ``examples/run_nue_12C_SF_1000.yml``, change the ``Cascade`` section:

.. code-block:: yaml

   Cascade:
     Run: True
     Interaction:
       Name: GeantInteractions
       GeantData: data/GeantData.hdf5
     Step: 0.04
     Probability: Cylinder

Setting ``Run: True`` activates the cascade. The key parameters are:

- **Interaction**: The interaction model used during propagation. ``GeantInteractions``
  uses tabulated nucleon-nucleon cross sections from Geant4 data.
- **Step**: The propagation step size in fm. Smaller values give higher accuracy.
- **Probability**: The model for computing interaction probability at each step.

For more sophisticated cascade configurations using the virtual-resonance model,
see the Cascade section of the :ref:`Run Card Structure` documentation, which
describes the ``NucleonNucleon``, ``DeltaInteraction``, and ``PionInteraction`` models.

What changes in the output
==========================

With the cascade enabled, events will contain additional vertices and particles
beyond the primary interaction:

- **Rescattered nucleons**: Nucleons that interacted with the nuclear medium during
  propagation, visible as cascade vertices.
- **Produced pions**: Secondary pions created through nucleon-nucleon or
  Delta-resonance interactions.
- **Absorbed particles**: Particles that were captured by the nucleus and do not
  appear in the final state.

The vertex status codes in the NuHepMC output distinguish primary interaction vertices
from cascade vertices, making it straightforward to separate the hard scattering from
final-state interactions in analysis.

Running with the cascade increases computation time per event because each hadronic
particle must be propagated step-by-step through the nucleus.
