.. _Physics Guide:

#############
Physics Guide
#############

This chapter explains what Achilles actually computes at each stage of an event,
so that users can make informed choices about run-card options. For the full
physics derivations and numerical validations, please refer to the publications
cited in :ref:`Introduction to Achilles`.

The Physics Guide is divided into the following sections:

- :doc:`physics/workflow` — how each event is generated end-to-end, and the
  integer codes for each nuclear interaction mode.
- :doc:`physics/nuclear-models` — the physics of each built-in nuclear model
  (Coherent, QESpectral, HyperonSpectral, FortranModel) and the three built-in
  Fortran models.
- :doc:`physics/form-factors` — the eight form factor categories, all available
  parametrisations, and how to change them.
- :doc:`physics/phase-space` — hadronic phase-space generators, Vegas and
  MultiChannel adaptive integration, and the Percentile unweighting algorithm.
- :doc:`physics/nucleus` — the physical meaning of each nucleus configuration
  field: binding energy, Fermi momentum, QMC configurations, and nuclear
  potential models.
- :doc:`physics/cascade` — intranuclear cascade algorithms, interaction
  probability models, in-medium corrections, and interaction models.
- :doc:`physics/ward-identity` — the Ward–Takahashi identity, why it is
  violated in the impulse approximation, and the four available gauge
  correction prescriptions.
