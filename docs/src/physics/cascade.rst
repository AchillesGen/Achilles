.. _cascade-physics:

***************************
Intranuclear Cascade Physics
***************************

After the hard interaction the struck nucleon — and any mesons or resonances
produced — must traverse the nuclear medium before being detected. The
intranuclear cascade (INC) models this propagation.

Step-based propagation
======================

The ``Base`` algorithm advances each active particle through the nucleus in
discrete steps of size ``Step`` (in fm). At each step the algorithm:

1. Computes the local nuclear density :math:`\rho(\boldsymbol{r})` at the
   particle's position.
2. Calculates the total nucleon–nucleon (or pion–nucleon) cross section
   :math:`\sigma` from the active interaction model.
3. Evaluates an interaction probability :math:`P_{\rm int}` from :math:`\sigma`
   and :math:`\rho` using the chosen ``Probability`` model.
4. Draws a random number to decide whether an interaction occurs.
5. If an interaction occurs, selects a target nucleon and applies the
   interaction model to determine the final state.
6. Applies Pauli blocking: forbids scattering into momentum states below the
   local Fermi momentum.

Particles that exit the nuclear volume escape; particles whose momentum falls
below the Fermi momentum are recaptured by the residual nucleus.

The ``MFP`` algorithm is intended to replicate the cascade algorithm in other event
generators: the path length to the next interaction would be drawn
from :math:`\ell \sim \text{Exp}(\lambda_{\rm mfp}^{-1})` where
:math:`\lambda_{\rm mfp} = 1/(\rho\,\sigma)`, which is more efficient when the
mean free path is large compared to the step size.

.. warning::

   The ``MFP`` algorithm is **not currently functional** — its interaction
   sampling routine is unimplemented and throws at runtime. Use the ``Base``
   algorithm for production runs.

Interaction probability models
================================

The ``Probability`` setting controls how the geometric cross section is converted
into an interaction probability at each step.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Model
     - Description
   * - ``Cylinder``
     - Each hadron is modelled as a hard cylinder of cross-sectional area
       :math:`\sigma`. An interaction occurs if the impact parameter between
       the propagating particle and a spectator falls within
       :math:`b < \sqrt{\sigma/\pi}`.
   * - ``Gaussian``
     - Each hadron is modelled with a Gaussian spatial distribution of width
       related to :math:`\sigma`. The interaction probability is the overlap
       integral of two Gaussians, smoothing the sharp cylinder edge.

In-medium corrections
======================

Nuclear matter modifies the effective NN cross section relative to its free-space
value. The ``InMedium`` setting controls the level of approximation.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Setting
     - Description
   * - ``None``
     - Free-space cross sections are used throughout. No density or Fermi-momentum
       dependence. Fastest option; adequate when cascade effects are a small
       correction.
   * - ``NonRelativistic``
     - The NN cross section is multiplied by a non-relativistic Pauli suppression
       factor that accounts for the fraction of final-state phase space blocked by
       the Fermi sea.
   * - ``Relativistic``
     - Reserved for a full relativistic in-medium modification (effective masses
       and mean-field potentials). **Not yet implemented:** selecting it is
       currently equivalent to ``None`` — the in-medium correction factor is
       applied only for ``NonRelativistic``.

Potential propagation
=====================

When ``PotentialProp: true``, particle trajectories are integrated with a
symplectic (Hamiltonian-preserving) integrator in the nuclear optical potential.
The potential bends particle paths and shifts their energy, allowing nucleons near
threshold to be captured into bound states. This is most relevant for low-energy
(sub-Fermi-momentum) nucleons. At higher energies the effect is small and
``PotentialProp: false`` is a good approximation.

Interaction models
==================

**NucleonNucleon**

Handles elastic and inelastic NN scattering. The ``GiBUU`` parametrisation uses
energy-dependent cross sections from the GiBUU transport model.

``ResonanceMode`` controls what happens when a :math:`\Delta` resonance is
produced in NN scattering:

- ``Decay`` — the :math:`\Delta` decays to :math:`N\pi` before the next cascade
  step. This is the simpler, faster option.
- ``Propagate`` — the :math:`\Delta` is propagated as an explicit particle
  through the nucleus and interacts via the ``DeltaInteraction`` model before
  decaying. Use ``Propagate`` together with ``DeltaInteraction`` in the
  interactions list.

``ExpSup`` is an exponential suppression factor applied at short inter-nucleon
distances to avoid unphysical hard scattering. Set to ``0.0`` to disable.

**DeltaInteraction**

Required when ``NucleonNucleon → ResonanceMode: Propagate``. Models
:math:`\Delta N` scattering and :math:`\Delta` propagation through nuclear matter.
``SWaveAbsorption: true`` enables s-wave :math:`\Delta` absorption
(:math:`\Delta N \to NN`).

**PionInteraction**

Models pion rescattering and absorption. Has two sub-models:

- ``HardScatter`` — quasi-elastic :math:`\pi N \to \pi N` scattering via
  ``MesonBaryonInteraction``.
- ``Absorption`` — pion absorption on a nucleon pair
  (:math:`\pi NN \to NN`) via ``PionAbsorptionOneStep``.

**ConstantInteraction**

Assigns user-specified constant cross sections to any two-body channel not
handled by the other models. This is used to give zero cross sections to photon
and hyperon channels so that they propagate freely through the nucleus.
