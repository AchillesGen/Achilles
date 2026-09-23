.. _form-factors-physics:

************
Form Factors
************

The ``FormFactorFile`` supplies one parametrisation for each of the eight form
factor *categories*. All categories are evaluated at every phase-space point, and
the appropriate subset is forwarded to the active nuclear model.

Categories and their physical role
====================================

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Category key
     - Computed quantities
     - Physical role
   * - ``vector``
     - :math:`F_{1p}`, :math:`F_{1n}`, :math:`F_{2p}`, :math:`F_{2n}`
     - Nucleon electromagnetic (Dirac and Pauli) form factors. Linear combinations
       :math:`F_1 = g_V^p F_{1p} + g_V^n F_{1n}` and similarly for :math:`F_2`
       are formed using the process-specific electroweak couplings.
   * - ``axial``
     - :math:`F_A`, :math:`F_{AP}`
     - Nucleon axial-vector form factor :math:`F_A(Q^2)` and the induced
       pseudoscalar :math:`F_{AP}(Q^2)`. Determines the parity-violating part
       of the weak-current cross section.
   * - ``coherent``
     - :math:`F_{\rm coh}`
     - Nuclear charge (coherent) form factor entering coherent scattering.
   * - ``resonancevector``
     - :math:`F_{\rm resV}`
     - Vector transition form factor for baryon resonance (e.g. :math:`\Delta`)
       excitation. Used by Fortran resonance models.
   * - ``resonanceaxial``
     - :math:`F_{\rm resA}`
     - Axial transition form factor for baryon resonance excitation.
   * - ``mecvector``
     - :math:`F_{\rm mecV3}`, :math:`F_{\rm mecV4}`, :math:`F_{\rm mecV5}`
     - Meson-exchange current vector contact terms. Parametrised as contact
       interactions with vector-meson propagator suppression.
   * - ``mecaxial``
     - :math:`F_{\rm mecA5}`
     - Meson-exchange current axial contact term involving the
       :math:`\Delta`-resonance.
   * - ``hyperon``
     - :math:`F_{1\Lambda}`, :math:`F_{2\Lambda}`, :math:`F_{A\Lambda}`,
       :math:`F_{1\Sigma^-}`, :math:`F_{2\Sigma^-}`, :math:`F_{A\Sigma^-}`,
       :math:`F_{1\Sigma^0}`, :math:`F_{2\Sigma^0}`, :math:`F_{A\Sigma^0}`
     - Nucleon-to-hyperon transition form factors for the three CC hyperon
       channels. Currently implemented as a placeholder (``Hyperon`` dummy).

Available parametrisations
============================

The table below lists every built-in parametrisation, the category it belongs to,
and a brief description.

.. list-table::
   :header-rows: 1
   :widths: 25 18 57

   * - Name
     - Category
     - Description
   * - ``Kelly``
     - ``vector``
     - Rational parametrisation of the Sachs form factors
       :math:`G_{Ep}`, :math:`G_{En}`, :math:`G_{Mp}`, :math:`G_{Mn}` from
       Kelly (2004). Recommended for quasielastic cross sections.
   * - ``BBBA``
     - ``vector``
     - Bradford–Budd–Bodek–Arrington parametrisation (2006) using rational
       polynomial fits to the Sachs form factors.
   * - ``ArringtonHill``
     - ``vector``
     - :math:`z`-expansion parametrisation of the Sachs form factors by
       Arrington and Hill (2006). Uses separate :math:`z`-expansion series for
       each Sachs form factor with a shared conformal parameter.
   * - ``VectorDipole``
     - ``vector``
     - Simple dipole approximation:
       :math:`G_D(Q^2) = (1 + Q^2/\lambda^2)^{-2}`.
       Useful for quick estimates; less accurate at large :math:`Q^2`.
   * - ``VectorDummy``
     - ``vector``
     - Sets :math:`F_{1p}=1`, all others zero. For debugging only.
   * - ``AxialZExpansion``
     - ``axial``
     - :math:`z`-expansion of the axial form factor :math:`F_A(Q^2)` and the
       strange axial contribution :math:`F_{As}(Q^2)`. The conformal variable
       :math:`z` maps the physical region onto the unit disk. Recommended
       parametrisation for the axial form factor.
   * - ``AxialDipole``
     - ``axial``
     - Dipole approximation
       :math:`F_A(Q^2) = g_A/(1+Q^2/M_A^2)^2`.
       The axial mass :math:`M_A` and the axial coupling :math:`g_A` are
       parameters.
   * - ``AxialDummy``
     - ``axial``
     - Sets :math:`F_A = F_{As} = 0`. For debugging only.
   * - ``Helm``
     - ``coherent``
     - Helm model of the nuclear charge form factor:
       :math:`F_{\rm coh}(Q^2) = 3 j_1(q R_0)/(q R_0)\,e^{-q^2 s^2/2}`.
       Suitable for medium and heavy nuclei. Parameters: skin thickness
       :math:`s` and nuclear radius derived from mass number :math:`A`.
   * - ``Lovato``
     - ``coherent``
     - *Ab initio* nuclear charge form factor from Lovato *et al.* for light
       nuclei (:math:`^{12}{\rm C}`). Uses a polynomial expansion fit to
       quantum Monte Carlo results.
   * - ``KN``
     - ``coherent``
     - Klein–Nystrand parametrisation.
   * - ``ResonanceVectorDummy``
     - ``resonancevector``
     - Constant :math:`F_{\rm resV} = {\tt resV}`. Placeholder until a full
       resonance form factor is implemented.
   * - ``ResonanceAxialDummy``
     - ``resonanceaxial``
     - Constant :math:`F_{\rm resA} = {\tt resA}`. Placeholder.
   * - ``MesonExchangeVector``
     - ``mecvector``
     - Serot–Walecka meson-exchange current vector contact terms with
       vector-meson propagator suppression
       :math:`\propto M_V^2/(M_V^2 + Q^2)`. Parameters: :math:`M_V^2` and
       three normalisation constants.
   * - ``MesonExchangeAxial``
     - ``mecaxial``
     - MEC axial contact term mediated by the :math:`\Delta` resonance.
       Parameters: :math:`M_{\Delta A}^2` and the :math:`C_5^A` normalisation.
   * - ``Hyperon``
     - ``hyperon``
     - Placeholder implementation. Returns zeros for all hyperon form factors
       until a proper parametrisation is added.

Changing a parametrisation
============================

To swap, for example, the vector form factor from ``Kelly`` to ``BBBA``, edit
the ``FormFactors.yml`` file:

.. code-block:: yaml

   # Change this line
   vector: BBBA

   # And add (or update) the BBBA parameter block
   BBBA:
     Mu Proton:  2.79278
     Mu Neutron: -1.91315
     NumeratorEp Params:   [1.0, -0.0578, 0.0, 0.0]
     DenominatorEp Params: [11.1, 13.6, 33.0, 0.0]
     ...

All other categories remain unchanged. Form factor variations are a common source
of theory uncertainty and should be explored when comparing to data.
