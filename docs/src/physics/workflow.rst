.. _event-workflow:

**************************
Event Generation Workflow
**************************

Each Achilles event passes through the stages shown below. The run-card sections
that control each stage are indicated in parentheses.

.. code-block:: text

   Beam sampling           (Beams)
         │
         ▼
   Process + mode selection (Processes, NuclearModels)
         │
         ▼
   Initial-state phase-space sampling
     – spectral function or Fermi gas drawn  (NuclearModels, Nuclei)
     – nucleon momentum p and removal energy E_r assigned
         │
         ▼
   Hard-scattering current calculation
     – hadronic transition current J^μ computed  (NuclearModels)
     – form factors evaluated at Q²              (FormFactorFile)
     – Ward gauge correction applied             (Ward)
         │
         ▼
   Cross-section weight assembled
     – leptonic current × hadronic current contracted
     – initial-state spectral-function weight multiplied in
         │
         ▼
   Intranuclear cascade  (Cascade)
     – struck nucleon and any produced hadrons
       stepped through the nuclear medium
         │
         ▼
   Event unweighting      (Options → Unweighting)
         │
         ▼
   Output                 (Main → Output)

**Beam sampling.** For each phase-space point the beam energy is drawn from the
configured flux. A ``Monochromatic`` beam always returns the same energy; a
``Spectrum`` beam samples from a histogram using the inverse-CDF method; a
``FlatFlux`` beam samples uniformly in energy.

**Process and mode selection.** The ``Processes`` list defines the leptonic
initial and final states. Combined with each ``NuclearModel``, Achilles enumerates
all allowed hadronic sub-processes (e.g. :math:`p \to p` for NC-QE,
:math:`p \to \Lambda` for CC-Hyperon). Each sub-process is integrated
independently and the results are summed incoherently.

**Initial-state phase-space sampling.** A nucleon is drawn from the nucleus
according to the nuclear model's phase-space generator
(``OneBodySpectral``, ``Coherent``, etc.). For spectral-function models this
samples a three-momentum :math:`\boldsymbol{p}` and removal energy :math:`E_r`
from the spectral function :math:`P(|\boldsymbol{p}|, E_r)`. The initial-state
weight for a nucleus with :math:`N_{p/n}` protons/neutrons is

.. math::

   w_{\rm IS} = N_{p/n}\,P(|\boldsymbol{p}|,\,E_r).

**Hard-scattering current.** The nuclear model computes the hadronic transition
current :math:`J^\mu` for the chosen hadronic subprocess, given the initial and
final nucleon momenta and the four-momentum transfer :math:`q^\mu`. Form factors
evaluated at :math:`Q^2 = -q^2` are passed into the current. The leptonic current
:math:`L_{\mu}` is assembled from the lepton kinematics (handled by the backend),
and the event weight receives a factor :math:`|L_{\mu} W^{\mu}|^2` where
:math:`W^{\mu}` is the hadronic current.

**Cascade.** If ``Cascade → Run: true``, struck nucleons and produced pions or
resonances are propagated through the nuclear medium. Particles that re-scatter,
get absorbed, or escape the nucleus are recorded accordingly.

**Unweighting.** After the integration phase converges (controlled by
``Options → Initialize → Accuracy``), Achilles generates unweighted events by
accept/reject on the event weight. The maximum weight used for rejection is
estimated from the ``Percentile`` of the weight distribution seen during
integration.


.. _nuclear-modes:

*************
Nuclear Modes
*************

Every nuclear model reports a mode integer that identifies the class of nuclear
interaction. This integer appears in the NuHepMC output and is used internally
to avoid assigning two models the same mode. The available modes are:

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Mode name
     - Integer
     - Description
   * - ``Coherent``
     - 1
     - The entire nucleus participates coherently. The nucleus recoils as a whole
       without break-up. Only neutral-current processes are currently allowed
       (leptonic charge = 0).
   * - ``Quasielastic``
     - 2
     - A single nucleon is struck and ejected. The residual nucleus is treated as
       a spectator. Allowed hadronic transitions:
       :math:`n \to n` and :math:`p \to p` (NC),
       :math:`n \to p` (CC :math:`\nu`) or :math:`p \to n` (CC :math:`\bar\nu`).
   * - ``MesonExchangeCurrent``
     - 3
     - Two correlated nucleons are struck simultaneously via a meson-exchange
       current. The hadronic state has two nucleons in both initial and final
       states (2p2h process). 16 spin states are summed. *Currently being implemented*
   * - ``Resonance``
     - 4
     - Single-pion production via baryon resonances (dominantly the
       :math:`\Delta(1232)`). The hadronic final state contains a nucleon and a
       pion, e.g. :math:`n \to p\,\pi^0`.
   * - ``ShallowInelastic``
     - 5
     - Shallow inelastic scattering. *Not yet implemented.*
   * - ``DeepInelastic``
     - 6
     - Deep inelastic scattering. *Not yet implemented.*
   * - ``Interference_QE_MEC``
     - 7
     - QE–MEC interference term. The hadronic state has one active nucleon and
       one spectator nucleon. Used together with the ``Quasielastic`` and
       ``MesonExchangeCurrent`` modes to capture cross-terms.
   * - ``Hyperon``
     - 8
     - Strange-quark hyperon production via charged-current interactions. Allowed
       transitions (charge +1):
       :math:`p \to \Lambda^0`, :math:`p \to \Sigma^0`, :math:`n \to \Sigma^-`.
