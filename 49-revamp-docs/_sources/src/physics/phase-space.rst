.. _phase-space-unweighting:

*************************************
Phase-Space Sampling and Unweighting
*************************************

Achilles uses adaptive Monte Carlo integration to evaluate the cross section and
generate events. This section describes the integration machinery, the
hadronic phase-space generators, and the unweighting algorithm.

Hadronic phase-space generators
================================

Each nuclear model declares the name of the phase-space generator it requires
via its ``PSName()`` method. The generator is responsible for mapping random
numbers to hadronic initial-state four-momenta. The three built-in generators
are:

.. list-table::
   :header-rows: 1
   :widths: 22 10 68

   * - Generator name
     - Dimensions
     - What it samples
   * - ``OneBodySpectral``
     - 4
     - Used by ``QESpectral``, ``HyperonSpectral``, and the Fortran QE and RES
       models. Draws one nucleon from the spectral function: nucleon momentum
       magnitude :math:`|\boldsymbol{p}|`, polar angle :math:`\theta`, azimuthal
       angle :math:`\phi`, and removal energy :math:`E_r`. Kinematic limits are
       set dynamically from the beam energy and nucleon mass.
   * - ``Coherent``
     - 0
     - Used by the ``Coherent`` model. The hadronic state is the nucleus at rest,
       so no random variables are needed.
   * - ``IntfSpectral``
     - 7
     - Used by the Fortran ``Intf_Spectral_Func`` model. Draws two nucleons: the
       active nucleon (3 variables: :math:`|\boldsymbol{p}_1|`, :math:`\theta_1`,
       :math:`\phi_1`) plus its removal energy :math:`E_{r,1}`, and the spectator
       nucleon (3 variables: :math:`|\boldsymbol{p}_2|`, :math:`\theta_2`,
       :math:`\phi_2` sampled uniformly up to 400 MeV).

Vegas adaptive integration
===========================

The innermost integration engine is **Vegas** :cite:`Lepage:1977sw`, which
adapts a multi-dimensional grid to importance-sample the integrand efficiently.

The integration proceeds in iterations. After each iteration the grid boundaries
are shifted to allocate more sampling density to bins that contributed the most
variance. This is controlled by the adaptation parameter ``alpha`` (default 1.5;
higher values produce more aggressive re-binning). The grid is periodically
refined by halving bin widths and doubling the number of calls per iteration.

Convergence is declared when both of the following are satisfied
simultaneously:

.. math::

   \frac{\sigma_I}{|\bar{I}|} < r_{\rm tol}
   \quad \text{and} \quad
   \sigma_I < a_{\rm tol}

where :math:`\bar{I}` is the current integral estimate, :math:`\sigma_I` is its
standard error, and :math:`r_{\rm tol}`, :math:`a_{\rm tol}` are the relative
and absolute tolerances. The ``Accuracy`` run-card parameter sets
:math:`r_{\rm tol}` (and :math:`a_{\rm tol}` to the same value).

MultiChannel integration
=========================

Achilles uses **MultiChannel** integration on top of Vegas. Several importance-
sampling channels (corresponding to different phase-space topologies or Feynman
diagrams) run in parallel, each with its own Vegas grid. After each iteration the
channel weights are updated: channels that explain more of the integrand's
variance receive larger weights, while channels with negligible contribution are
suppressed (but never below a minimum floor ``min_alpha = 1e-5`` to prevent
permanent elimination).

The combined estimator is:

.. math::

   I = \sum_i \alpha_i \, I_i, \qquad \alpha_i \geq 0,\;\sum_i\alpha_i = 1

where :math:`\alpha_i` are the adaptive channel weights and :math:`I_i` is the
integral from channel :math:`i`. The weights evolve automatically and do not
require user input. The ``Accuracy`` parameter controls the termination of the
multichannel optimisation phase via the same relative-tolerance criterion as
Vegas.

Options settings and their effect
===================================

The ``Options`` run-card section controls both the random seed and the integration
accuracy. For details of the YAML syntax see the run-card documentation; the
physics meaning of each setting is:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Setting
     - Effect
   * - ``Initialize → Seed``
     - Seeds the global pseudo-random number generator before any integration
       or sampling begins. Using the same seed with the same run card exactly
       reproduces the sequence of random numbers and therefore the same sequence
       of accepted events (up to floating-point portability). Change the seed to
       generate statistically independent runs.
   * - ``Initialize → Accuracy``
     - Sets the relative (and absolute) tolerance :math:`r_{\rm tol}` for the
       convergence test. Smaller values require more Vegas iterations and function
       evaluations before the integration is considered converged, but produce a
       more reliable estimate of the total cross section. Typical values:
       ``1e-2`` (fast, for exploratory runs), ``1e-3`` (production quality).

Unweighted events and the Percentile algorithm
================================================

After the integration phase, Achilles generates *unweighted* events — events
that each represent the same "amount" of physics and can be treated as equally
likely. This is done by accept/reject on the event weight :math:`w`:

.. math::

   P(\text{accept}) = \frac{|w|}{w_{\rm max}}

where :math:`w_{\rm max}` is a chosen maximum weight. Events with
:math:`|w| > w_{\rm max}` are accepted with probability 1 but assigned a
residual weight :math:`|w|/w_{\rm max}` (they remain partially weighted).

The ``Percentile`` algorithm sets :math:`w_{\rm max}` to the :math:`p`-th
percentile of all event weights seen during integration:

- **High percentile** (e.g. 99): :math:`w_{\rm max}` is large, so almost all
  events are accepted, but rare large-weight events may survive as partially
  weighted outliers.
- **Low percentile** (e.g. 80): :math:`w_{\rm max}` is small, so more events
  are rejected outright, giving a more uniform sample at the cost of generation
  efficiency.

The percentile is tracked exactly using an online algorithm as integration
proceeds, so no second pass over the data is needed.

The ``percentile`` run-card parameter (0–100) directly sets :math:`p`.
A value of ``99`` is recommended for most production runs: it discards the
top 1 % of weights, which are typically integration-grid artefacts, while
keeping 99 % of events fully unweighted.

The unweighting efficiency (fraction of phase-space points that produce an
accepted event) is printed in the Achilles log after event generation.
