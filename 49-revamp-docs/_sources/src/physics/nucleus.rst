.. _nucleus-config-physics:

**********************
Nucleus Configuration
**********************

The ``Nuclei`` run-card section defines the nuclear target. This section
explains the physical meaning of each field in a nucleus file.

Binding energy and Fermi momentum
===================================

``Binding`` is the average binding energy per nucleon in MeV. It is used as a
global offset when computing removal energies in models that do not read a
spectral function (e.g. a simple Fermi-gas approximation).

``Fermi Momentum`` is the Fermi momentum :math:`k_F` in MeV, used as the
upper limit of the occupied momentum states for Pauli blocking during the
cascade. For a local Fermi gas the effective :math:`k_F` varies with position
according to the local density; the ``Fermi Momentum`` value is used as a
global (average) fallback.

Nuclear density and configuration files
=========================================

The ``Density`` sub-node provides the spatial description of the nucleus.

``Function: configuration`` activates the quantum Monte Carlo (QMC)
configuration sampling mode. Instead of using a smooth analytic density, each
event draws a snapshot of nucleon positions from a pre-generated library of
nuclear ground-state configurations computed by QMC methods. This captures
short-range correlations and clustering effects that a smooth density misses.

The configurations are stored in a gzip-compressed text file given by the
``Configs`` key. Achilles reads one configuration per event (cycling through
the file), so the number of distinct configurations sets an upper limit on
the number of statistically independent initial nuclear states.

The ``ProtonFile`` and ``NeutronFile`` keys provide tabulated radial density
profiles :math:`\rho_{p/n}(r)`, used both to evaluate the local nuclear density
during the cascade and to draw the nuclear potential.

``FilePotential`` provides the real part of an optical potential (used when
``PotentialProp: true`` in the cascade).

Fermi-gas parameters
======================

The ``FermiGas`` sub-node controls corrections to the basic Fermi-gas picture
that can optionally be applied within the cascade.

``Type: Local`` uses a local Fermi gas, where the Fermi momentum varies with
the local nuclear density:

.. math::

   k_F(\boldsymbol{r}) = \left(\frac{3\pi^2}{2}\,\rho(\boldsymbol{r})\right)^{1/3}

``Correlated`` and the SRC parameters add short-range correlation (SRC)
corrections to the nucleon momentum distribution. When ``Correlated: true`` a
fraction ``SRCfraction`` of nucleons are assigned high-momentum tails extending
beyond :math:`k_F`, parametrised by an exponential with scale ``LambdaSRC``
(in GeV). This is important for nucleon momenta above ~300 MeV where the
spectral function departs significantly from a simple Fermi gas. The default
values ``SRCfraction: 0.2`` and ``LambdaSRC: 2.75`` are tuned to reproduce
electron-scattering data.

Nuclear potential models
=========================

The ``Potential`` sub-node selects the nuclear optical potential used for
cascade propagation (when ``PotentialProp: true``) and Pauli blocking.

**``Wiringa``**

A non-relativistic density-dependent potential :cite:`Wiringa:1988jt`. The
central potential is:

.. math::

   V(p_{\rm lab}, r) = \alpha(\rho) + \frac{\beta(\rho)}{1 + (p_{\rm lab}/\lambda(\rho))^2}

where :math:`\rho = \rho(r)/\rho_0` is the local density normalised to
saturation density :math:`\rho_0` (set by ``r0``, default 0.16 fm\ :sup:`-3`),
and the coefficients are:

.. math::

   \alpha(\rho) &= 15.52\,\rho + 24.93\,\rho^2 \;\text{MeV} \\
   \beta(\rho)  &= -116\,\rho \;\text{MeV} \\
   \lambda(\rho)&= (3.29 - 0.373\,\rho)\,\hbar c

The Hamiltonian used for propagation is non-relativistic:
:math:`H = m_N + p^2/(2m_N) + V`.

**``Cooper``**

A relativistic optical potential :cite:`Cooper:2009zza` in the relativistic
impulse approximation. It provides separate real and imaginary parts for the
vector and scalar self-energies (:math:`\Sigma_V`, :math:`\Sigma_S`) as
polynomial functions of the local density and the projectile energy. This
is the most comprehensive potential available and is recommended for heavy
nuclei (:math:`A \geq 12`).

**``Schroedinger``**

A Schrödinger-equivalent potential derived from the Cooper relativistic optical
model by reducing the Dirac equation to a second-order Schrödinger form with
Darwin correction terms. The ``Mode`` parameter selects the reduced-mass
prescription used in the reduction:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Mode
     - Reduced-mass prescription
   * - 2
     - :math:`\mu = E_{\rm cm,p}\,E_{\rm cm,t}/\sqrt{s}`  (kinematic)
   * - 3
     - :math:`\mu = m_p\,m_t/(m_p+m_t)`  (standard non-relativistic)
   * - 4
     - :math:`\mu = E_{\rm cm,p}`  (projectile energy only)
   * - 5
     - :math:`\mu = m_p`  (projectile mass only)

Mode 3 is recommended for cascade propagation of nucleons in medium-mass nuclei.
The Hamiltonian is non-relativistic: :math:`H = m_N + p^2/(2m_N) + V_{\rm central}`.
