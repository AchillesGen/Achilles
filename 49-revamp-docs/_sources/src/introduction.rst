.. _Introduction:

############
Introduction
############

Achilles (**A** **CHI**\ cag\o **L**\ and **L**\ epton **E**\ vent **S**\ imulator) is a 
Monte Carlo event generator for the simulation of lepton-nucleus and neutrino-nucleus interactions.
This manual provides information to help users understand and apply Achilles for their physics studies.

The general organization of the manual is as follows. First, the generator and the overarching goals
are introduced. Secondly, a detailed description of installing and running the program are outlined.
Finally, methods of developing plugins to extend the functionality are covered, along with detailed
documentation of the needed interfaces.
The various options and parameters specifying how the program is to be compiled and run, and their
meanings are explained.

Note: This document is meant as a tool for using and running the code, therefore,
a complete detailed description of the physics content of Achilles is not available here.
For those interested in the physics details, the authors refer the reader to the following publications,
:cite:`Isaacson:2025cnk`, :cite:`Isaacson:2022cwh`, :cite:`Isaacson:2021xty`, and :cite:`Isaacson:2020wlx`.

.. contents::
   :local:

.. _Introduction to Achilles:

************************
Introduction to Achilles
************************

Achilles is a Monte Carlo event generator that provides fully exclusive final states in simulations for
lepton-nucleus and neutrino-nucleus interactions.
The produced events may be passed into detector simulations used by the various experiments, or analyzed
for comparisons to data.
The Achilles code consists of a main C++ code, with extensions to include nuclear models implemented in modern
Fortran90 (see :ref:`Fortran Interface <fortran-interface>` for more details).

The main goals of the Achilles event generator are:

1. Theory driven: ensure that the guiding principles are to reproduce our understanding of the
   Standard Model with appropriate uncertainties.
2. Leverage experiences from the LHC event generator community
3. Develop a modular and extensible neutrino event generator
4. Provide automated Beyond the Standard Model calculations for neutrino experiments
5. Fully evaluate theory uncertainties

In addition to lepton-nucleus and neutrino-nucleus run modes, Achilles provides a setup to benchmark the
intranuclear cascade. This enables running nucleon-nucleus and pion-nucleus reaction cross-sections as well.
The list of allowed nuclear physics interactions currently includes coherent scattering, quasielastic scattering,
and resonance production. Ongoing efforts are aimed at developing models for coherent single pion production,
meson exchange currents (MEC), shallow inelastic scattering (SIS), and deep inelastic scattering (DIS).
In terms of allowed leptonic initial and final states, it is possible to simulate any set of colorless particles,
including those beyond the Standard Model (SM) through the use of an :ref:`interface <sherpa-interface>` to
`Sherpa`_, the UFO output format :cite:`Degrande:2011ua,Darme:2023jdn` from
`FeynRules`_:cite:`Christensen:2008py,Christensen:2009jx,Alloul:2013bka`.

This manual aims to provide all the needed information to get started with Achilles as quickly as possible.
It lists all available options for controlling the physics aspects of the simulation.
It does not describe the details of the physics simulated, and only minor details of the underlying structure
as it pertains to enabling user extensions.

The `MCnet Guidelines
<https://www.montecarlonet.org/publications_guidelines/>`_
apply to Achilles. Therefore, you are kindly reminded that the guidelines recommend the citation of all
works used for their continued development. We kindly ask that you cite :cite:`Isaacson:2025cnk,Isaacson:2022cwh`
for all works using Achilles. Additionally, the Achilles generator outputs a bibtex file containing additional
works that should be cited based upon the run configuration used to generator your events, and also included
in the `NuHepMC`_:cite:`Gardiner:2023ejq` output format.

The Achilles authors strongly recommend a detailed review of the manuals and many excellent publications
on different aspects of event generation within the neutrino community and beyond. An excellent resource
for those interested in seeing the current state and future goals of event generators in the community
is the white paper written for Snowmass :cite:`Campbell:2022qmc`.

The remainder of this manual is organized as follows: in :ref:`Basic Structure` the modular structure
is highlighted. :ref:`Getting Started` contains information and instructions on obtaining, building,
installing, and running the event generator. The :ref:`Run Configuration` is discussed, and the ability
to control the parameters and physics within the generator are euclidated. Finally, the :ref:`Interfaces`
for extending the code or to external tools are explained.

**Note**: The construction of Monte Carlo event generators requires several assumptions, approximations,
and simplifications. It is therefore imperative that the results of event generators should be verified and
cross-checked with results from other programs, and they should be interpreted with care and common sense.
For a list of other neutrino event generators, please see the following sub-section.


.. _Other Generators:

----------------
Other Generators
----------------

Other neutrino event generators that exist with various states of support, extensibility, goals, and level of open-access include:

* `GENIE <https://genie-mc.github.io/>`_
* `GiBUU <https://gibuu.hepforge.org/>`_
* `Marley <https://www.marleygen.org/>`_
* `NuWro <https://nuwro.github.io/user-guide/>`_
* NEUT (Currently closed source)

.. _Sherpa: https://sherpa-team.gitlab.io/
.. _Feynrules: http://feynrules.irmp.ucl.ac.be/
.. _NuHepMC: https://github.com/nuhepMC/Spec
