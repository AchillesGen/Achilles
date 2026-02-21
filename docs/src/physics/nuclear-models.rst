.. _nuclear-models-physics:

**************
Nuclear Models
**************

.. contents::
   :local:

``Coherent``
============

In coherent scattering the whole nucleus participates and the nucleus remains in
its ground state. The hadronic current is proportional to the sum of initial and
final nuclear four-momenta weighted by the nuclear form factor:

.. math::

   J^\mu_{\rm coh} = \bigl(p_A^\mu + p_A^{\prime\,\mu}\bigr)\,F_{\rm coh}(Q^2)

where :math:`F_{\rm coh}(Q^2)` is the coherent (nuclear charge) form factor. Only
the ``coherent`` entry in the form factor file is used. The initial-state weight is
unity (the nucleus is always in its ground state).

``QESpectral``
==============

Quasi-elastic scattering is treated in the impulse approximation: a single
off-shell nucleon with three-momentum :math:`\boldsymbol{p}` and removal energy
:math:`E_r = m_N - E` is struck by the virtual boson and ejected from the nucleus.

**Initial-state sampling.** The joint probability for finding a nucleon with
momentum :math:`|\boldsymbol{p}|` and removal energy :math:`E_r` is given by the
nucleon spectral function:

.. math::

   w_{\rm IS} = N_{p/n}\,P(|\boldsymbol{p}|,\,E_r)

where :math:`N_{p/n}` is the number of protons or neutrons and :math:`P` is read
from ``SpectralP`` / ``SpectralN``. The spectral function is stored as a tabulated
2D array in momentum and removal energy and is interpolated at run time.

**Hadronic current.** The current uses the standard nucleon vertex with Dirac
(:math:`F_1`), Pauli (:math:`F_2`), and axial (:math:`F_A`, :math:`F_{AP}`)
form factors:

.. math::

   J^\mu = \bar{u}(p^\prime)\Bigl[
       F_1\,\gamma^\mu
       + i\,\frac{F_2}{2m_N}\,\sigma^{\mu\nu}q_\nu
       + F_A\,\gamma^\mu\gamma^5
       + \frac{F_{AP}}{m_N}\,\gamma^5 q^\mu
   \Bigr]u(p)

where :math:`u(p)` and :math:`\bar{u}(p^\prime)` are free Dirac spinors for the
initial and final nucleons. The initial nucleon energy is set to the free on-shell
value :math:`E_{\rm free} = \sqrt{|\boldsymbol{p}|^2 + m_N^2}`, while the energy
transfer is shifted by the off-shell correction
:math:`\omega \to \omega + E - E_{\rm free}`.

All four spin-state combinations :math:`(s, s^\prime) \in \{-1,+1\}^2` are summed
incoherently, giving 4 entries in the current array.

``HyperonSpectral``
===================

Hyperon production follows the same spectral-function formalism as ``QESpectral``
for the initial state. The hadronic current uses separate Dirac, Pauli, and axial
form factors for each hyperon channel (:math:`\Lambda^0`, :math:`\Sigma^0`,
:math:`\Sigma^-`), taken from the ``hyperon`` block of the form factor file.
The outgoing baryon mass in the Pauli term is replaced by the hyperon mass.

``FortranModel``
================

``FortranModel`` acts as a bridge between the C++ framework and a user-supplied
Fortran nuclear model. The C++ class calls the Fortran routines described in
:ref:`fortran-interface` to obtain the hadronic current and initial-state weight.
The nuclear mode and phase-space generator are also determined by the Fortran model
at runtime, so a ``FortranModel`` can implement any of the modes listed in
:ref:`nuclear-modes`.

Three Fortran models ship with Achilles and are registered automatically at startup.
They are summarised below. For all three the run-card ``Model`` key must be set to
``FortranModel`` and the ``Name`` key must match the model's registered name exactly.

.. list-table:: Built-in Fortran models at a glance
   :header-rows: 1
   :widths: 28 22 10 10 30

   * - Registered name (``Name``)
     - Source file
     - Mode
     - Frame
     - InspireHEP
   * - ``QE_Spectral_Func``
     - ``qe_spectral_model.f90``
     - Quasielastic (2)
     - Lab (0)
     - ``Rocco:2018mwt``
   * - ``RES_Spectral_Func``
     - ``res_spectral_model.f90``
     - Resonance (4)
     - QZ (1)
     - ``Nakamura:2015rta``, ``Rocco:2019gfb``
   * - ``Intf_Spectral_Func``
     - ``intf_spectral_model.f90``
     - Interference QE–MEC (7)
     - Lab (0)
     - ``Lovato:2023khk``

``QE_Spectral_Func`` — Quasielastic spectral function
------------------------------------------------------

This is the Fortran implementation of the quasi-elastic spectral-function model,
equivalent in physics to the C++ ``QESpectral`` model. It samples a nucleon from
the nuclear spectral function and computes the one-body hadronic current with
Dirac, Pauli, and axial form factors.

**Physics.** The hadronic current is:

.. math::

   J^\mu = \bar{u}(p^\prime)\Bigl[
       F_1\,\gamma^\mu
       + i\,\frac{F_2}{2m_N}\,\sigma^{\mu\nu}q_\nu
       + F_A\,\gamma^\mu\gamma^5
       + \frac{F_{AP}}{m_N}\,\gamma^5 q^\mu
   \Bigr]u(p)

with the off-shell energy correction :math:`\omega \to \omega + E - E_{\rm free}`
applied to the energy transfer. Four spin combinations are summed (``nspin = 4``).
The initial-state weight is :math:`w = N_{p/n}\,P(|\boldsymbol{p}|, E_r)`.

**Form factors used:** ``F1``, ``F2``, ``FA``, ``FAP`` (from the ``vector`` and
``axial`` categories of the form factor file).

**Configuration.** The ``ConfigFile`` is a plain-text file with two lines giving
the paths to the proton and neutron spectral function data files:

.. code-block:: text

   data/Spectral_Functions/pke12p_tot.data
   data/Spectral_Functions/pke12n_tot.data

Pre-supplied config files: ``data/info_C12_pke.data``, ``data/info_O16_pke.data``,
``data/info_Ar40_pke.data``.

**Run-card example:**

.. code-block:: yaml

   NuclearModels:
     - NuclearModel:
         Model: FortranModel
         Name: QE_Spectral_Func
         ConfigFile: data/info_C12_pke.data
         FormFactorFile: FormFactors.yml
         Ward: None

``RES_Spectral_Func`` — Resonance spectral function
-----------------------------------------------------

Single-pion production through baryon resonances (primarily the
:math:`\Delta(1232)`), combined with a spectral-function description of the
nuclear initial state. All pion–nucleon channels consistent with charge
conservation are included automatically.

**Physics.** The hadronic state has one incoming nucleon (with spectral-function
weight) and one outgoing nucleon plus one pion. The current is computed in the
frame where :math:`\boldsymbol{q}` is aligned along the :math:`z`-axis
(``frame = QZ``), which simplifies the angular integration. The vector coupling
is set by the ``FResV`` form factor; the axial coupling is set by ``FResA``. The
``has_axial`` flag is ``true`` whenever ``FResA`` is non-zero, enabling axial
contributions for weak interactions while switching them off for purely
electromagnetic processes.

**Form factors used:** ``FResV``, ``FResA`` (from the ``resonancevector`` and
``resonanceaxial`` categories).

**Configuration.** The ``ConfigFile`` uses the same two-line format as
``QE_Spectral_Func``, supplying proton and neutron spectral function paths:

.. code-block:: text

   data/Spectral_Functions/pke12p_tot.data
   data/Spectral_Functions/pke12n_tot.data

The same pre-supplied config files (``data/info_C12_pke.data``, etc.) can be
re-used. No ``ModelParamsFile`` is needed.

**Run-card example:**

.. code-block:: yaml

   NuclearModels:
     - NuclearModel:
         Model: FortranModel
         Name: RES_Spectral_Func
         ConfigFile: data/info_C12_pke.data
         FormFactorFile: FormFactors.yml
         Ward: None

``Intf_Spectral_Func`` — QE–MEC interference spectral function
---------------------------------------------------------------

This model computes the *interference* between the one-body quasi-elastic
amplitude and the two-body meson-exchange current (MEC) amplitude. It
corresponds to nuclear mode ``Interference_QE_MEC`` (7) and requires its own
phase-space generator (``IntfSpectral``), which samples both an active nucleon
and a spectator nucleon from the nucleus.

**Physics.** The nuclear current is the cross-term between the one-body current
:math:`J^\mu_{1b}` (computed as in ``QE_Spectral_Func``) and the two-body MEC
current :math:`J^\mu_{2b}`. The two-body current receives contributions from
both pion-exchange and :math:`\Delta`-exchange diagrams, in direct and exchange
topologies:

.. math::

   J^\mu_{2b} = \frac{1}{2 E_{p_2}}\Bigl(
       J^\mu_{\pi,{\rm dir}} + J^\mu_{\Delta,{\rm dir}}
       + J^\mu_{\pi,{\rm exc}} + J^\mu_{\Delta,{\rm exc}}
   \Bigr)

where :math:`E_{p_2}` is the energy of the spectator nucleon.

Because the QE and MEC amplitudes share the same hadronic phase-space point,
the model alternates between returning the one-body current and the two-body
current on successive calls (controlled by an internal flag ``compute_1body``).
This paired-call structure is handled transparently by the Achilles framework.

**Initial-state weight.** The weight is the product of two spectral-function
factors — one for the active (hole) nucleon and one for the spectator:

.. math::

   w = \mathcal{N}\,P(|\boldsymbol{p}_1|, E_{r,1}) \;\times\;
       \mathcal{N}\,P(|\boldsymbol{p}_2|)

where :math:`\mathcal{N}` is the spectral function normalisation and the
spectator factor is obtained by integrating over removal energy.

**Form factors used:** ``F1``, ``F2``, ``FA``, ``FAP`` (one-body part);
``FMecV3``, ``FMecV4``, ``FMecV5``, ``FMecA5``, ``FPiEM`` (two-body MEC part).

**Configuration.** The ``ConfigFile`` lists the *mean-field* proton and neutron
spectral function files (distinct from the total spectral functions used in the
QE model):

.. code-block:: text

   data/Spectral_Functions/pke12p_MF.data
   data/Spectral_Functions/pke12n_MF.data

Pre-supplied config files: ``data/info_intf_C12.data``, ``data/info_intf_Ar40.data``.

A ``ModelParamsFile`` is **required** for this model. It supplies the hadronic
coupling constants and momentum cut-offs for the MEC vertex functions:

.. code-block:: yaml

   Parameters:
     fpind:  0.54      # f_{pi N Delta} coupling
     fstar:  2.13      # f* Delta coupling
     fpinn2: 1.01789   # f_{pi N N}^2 coupling
     ga:     1.26      # axial coupling g_A
     lpi:    1300.0    # Lambda_pi cut-off (MeV)
     lpind:  1150.0    # Lambda_{pi N Delta} cut-off (MeV)

The default parameter file ``data/intf_params.yml`` provides the values above,
which were used in :cite:`Lovato:2023khk`.

**Run-card example:**

.. code-block:: yaml

   NuclearModels:
     - NuclearModel:
         Model: FortranModel
         Name: Intf_Spectral_Func
         ConfigFile: data/info_intf_C12.data
         FormFactorFile: FormFactors.yml
         ModelParamsFile: data/intf_params.yml
         Ward: None
