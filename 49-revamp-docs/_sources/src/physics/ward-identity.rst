.. _ward-identity:

***********************************
Ward Identity and Gauge Corrections
***********************************

The Ward–Takahashi identity (WTI) of QED/QFD requires the hadronic current to
be conserved:

.. math::

   q_\mu J^\mu = 0

where :math:`q^\mu` is the four-momentum transfer. For free on-shell nucleons
the WTI is satisfied automatically when the correct form factors are used.
However, the nuclear off-shell treatment in the impulse approximation introduces
violations: the nucleon inside the nucleus has a modified energy
(:math:`E \ne \sqrt{|\boldsymbol{p}|^2+m_N^2}`) but its vertex is still
evaluated as if it were on-shell. This mismatch means
:math:`q_\mu J^\mu \ne 0` in general.

The ``Ward`` key selects one of four prescriptions for handling this:

``None``
========

No correction is applied. The current returned by the nuclear model is passed
to the cross-section calculation as-is.

**When to use:** When the nuclear model already satisfies the WTI internally
(e.g. some fully relativistic models), or when you want to directly measure the
size of the WTI violation as a systematic uncertainty.

``Coulomb``
===========

The longitudinal spatial component of the current along :math:`\boldsymbol{q}`
is fixed by the time (charge) component via the continuity equation:

.. math::

   J_L \;\to\; \frac{\omega}{|\boldsymbol{q}|}\,J^0

where :math:`\omega = q^0` is the energy transfer and
:math:`J_L = \hat{q}\cdot\boldsymbol{J}` is the longitudinal component. The
transverse spatial components and :math:`J^0` are left unchanged.

**When to use:** When the time component of the current is well-constrained by
the nuclear model (as in non-relativistic spectral-function calculations) and
the spatial current needs to be made consistent with charge conservation. This
is the most commonly used correction for quasi-elastic spectral-function models.

``Weyl``
========

The time component is fixed by the longitudinal spatial component:

.. math::

   J^0 \;\to\; \frac{|\boldsymbol{q}|}{\omega}\,J_L

This is complementary to the Coulomb prescription. The spatial components are
left unchanged, and :math:`J^0` is derived from them.

**When to use:** When the spatial components of the current are considered more
reliable than the time component.

``Landau``
==========

The full covariant Lorenz (Landau) gauge condition is applied. The current is
projected transverse to :math:`q^\mu` by subtracting its longitudinal part:

.. math::

   J^\mu \;\to\; J^\mu - \frac{J \cdot q}{Q^2}\,q^\mu

where :math:`Q^2 = -q^2 > 0`. This is a manifestly covariant operation.

**When to use:** When a fully covariant treatment is required, for example in
relativistic nuclear models where Lorentz covariance must be maintained
explicitly.

.. note::

   The choice of Ward correction does not change the physics for a model that
   already satisfies the WTI; it only affects the result when the WTI is
   violated. Comparing ``Ward: None`` with the other options therefore provides
   an estimate of the systematic uncertainty from nuclear off-shell effects.
   The default ``Ward: None`` is appropriate for a first exploration; for
   quantitative comparisons to data ``Ward: Coulomb`` is recommended for
   spectral-function models.
