.. _Interfaces

##########
Interfaces
##########

.. contents::
   :local:

.. _fortran-interface:

Fortran Interface
=================

Achilles provides a Fortran interface that allows nuclear models written in Fortran to be
used directly within the C++ framework. The interface is structured in three layers:

1. An abstract Fortran base type (``model``) that a user implements.
2. A Fortran factory (``ModelFactory``) that manages registration and instantiation of models.
3. A set of ISO C-bound procedures (``nuclear_model_interface``) that bridge Fortran to the
   C++ ``FortranModel`` class.

The overall flow when running a simulation is:

- At startup, ``FortranModel::RegisterModels()`` calls the Fortran ``RegisterAll()`` routine,
  which registers every available Fortran model with the factory.
- When the YAML run card requests a Fortran-backed nuclear model, ``FortranModel`` calls
  ``CreateModel`` to instantiate it and ``InitModel`` to configure it from the supplied
  configuration file and parameter map.
- During event generation ``FortranModel`` delegates the current calculation and initial-state
  weighting to the Fortran model via ``GetCurrents`` and ``GetInitialStateWeight``.

The ``model`` Abstract Type
---------------------------

*Defined in* ``src/Achilles/fortran/nuclear_model.f90``

Every Fortran nuclear model must extend the abstract ``model`` type and implement all of
its deferred procedures.

.. code-block:: fortran

   type, abstract :: model
       contains
           procedure(nm_init),     deferred          :: init
           procedure(nm_cleanup),  deferred          :: cleanup
           procedure(nm_mode),     deferred          :: mode
           procedure                                 :: frame => nm_frame
           procedure(nm_name),     deferred, nopass  :: model_name
           procedure(nm_psname),   deferred          :: ps_name
           procedure(nm_currents), deferred          :: currents
           procedure(nm_init_wgt), deferred          :: init_wgt
           procedure(nm_inspirehep), deferred, nopass :: inspirehep
   end type model

The table below describes each deferred procedure together with its signature and purpose.

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Procedure
     - Signature
     - Description
   * - ``init``
     - ``logical function init(self, filename, params)``
       where ``filename`` is a Fortran string and ``params`` is a ``map`` of
       ``character``-keyed ``real(c_double)`` values.
     - Load model-specific data from ``filename`` and apply any numeric parameters
       from ``params``. Return ``.true.`` on success.
   * - ``cleanup``
     - ``subroutine cleanup(self)``
     - Deallocate any resources allocated during ``init``.
   * - ``mode``
     - ``integer function mode(self)``
     - Return an integer identifying the nuclear interaction mode (e.g. quasi-elastic,
       resonance). The value must match the ``NuclearMode`` enum in the C++ layer.
   * - ``frame``
     - ``integer function frame(self)``
     - Return an integer identifying the preferred reference frame. The default
       implementation returns ``0`` (lab frame). Override only when a different frame
       is required.
   * - ``model_name``
     - ``character(len=:), allocatable function model_name()`` (``nopass``)
     - Return the unique string name used to identify this model in the run card and
       factory. Must match the name passed to ``register_model``.
   * - ``ps_name``
     - ``character(len=:), allocatable function ps_name(self)``
     - Return the name of the phase-space generator required by this model (e.g.
       ``"QESpectral"``).
   * - ``currents``
     - ``subroutine currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout,``
       ``pids_spect, mom_spect, nspect, qvec, ff, cur, nspin, nlorentz)``
     - Compute the hadronic transition current. See :ref:`currents-signature` below for
       the full argument list.
   * - ``init_wgt``
     - ``real(c_double) function init_wgt(self, pids_in, mom_in, nin, pids_spect,``
       ``mom_spect, nspect, nproton, nneutron)``
     - Return the weight from the nuclear initial state (e.g. spectral function or Fermi
       gas distribution). See :ref:`init-wgt-signature` below.
   * - ``inspirehep``
     - ``character(len=:), allocatable function inspirehep()`` (``nopass``)
     - Return the InspireHEP record identifier (e.g. ``"1234567"``) for the paper
       describing this model. Used for citations in the output.

.. _currents-signature:

``currents`` Argument List
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   subroutine nm_currents(self, pids_in, mom_in, nin, &
                          pids_out, mom_out, nout, &
                          pids_spect, mom_spect, nspect, &
                          qvec, ff, cur, nspin, nlorentz)
       class(model),         intent(inout) :: self
       integer(c_size_t),    intent(in), value :: nin, nout, nspect, nspin, nlorentz
       integer(c_long),      dimension(nin),    intent(in) :: pids_in
       integer(c_long),      dimension(nout),   intent(in) :: pids_out
       integer(c_long),      dimension(nspect), intent(in) :: pids_spect
       type(fourvector),     dimension(nin),    intent(in) :: mom_in
       type(fourvector),     dimension(nout),   intent(in) :: mom_out
       type(fourvector),     dimension(nspect), intent(in) :: mom_spect
       type(fourvector),               intent(in)  :: qvec
       type(complex_map),              intent(in)  :: ff
       complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur

``ff`` is a map from form-factor name strings to complex values, pre-evaluated at the
momentum transfer :math:`Q^2`. ``cur`` is filled column-by-column: ``cur(mu, spin)`` is
the :math:`\mu`-th Lorentz component for spin state ``spin``.

.. _init-wgt-signature:

``init_wgt`` Argument List
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   function nm_init_wgt(self, pids_in, mom_in, nin, &
                        pids_spect, mom_spect, nspect, &
                        nproton, nneutron) result(wgt)
       class(model),       intent(inout) :: self
       integer(c_size_t),  intent(in), value :: nin, nspect, nproton, nneutron
       integer(c_long),    dimension(nin),    intent(in) :: pids_in
       integer(c_long),    dimension(nspect), intent(in) :: pids_spect
       type(fourvector),   dimension(nin),    intent(in) :: mom_in
       type(fourvector),   dimension(nspect), intent(in) :: mom_spect
       real(c_double) :: wgt

``nproton`` and ``nneutron`` give the proton and neutron number of the target nucleus.
The return value is the probability weight for the sampled initial-state configuration.

The ``ModelFactory`` Type
-------------------------

*Defined in* ``src/Achilles/fortran/nuclear_model_factory.f90``

``ModelFactory`` maintains a list of named model constructors and creates model instances
on demand.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Procedure
     - Description
   * - ``init(self)``
     - Initialise the factory (nullify the internal model pointer).
   * - ``final(self)``
     - Deallocate and nullify the stored model pointer.
   * - ``register_model(self, name, construct)``
     - Register a new model. ``name`` is the unique string identifier (up to 20
       characters); ``construct`` is a ``nopass`` function pointer with signature
       ``class(model), pointer function constructor()``. Halts if ``name`` is already
       registered.
   * - ``create_model(self, name)``
     - Return a pointer to a newly-created model matching ``name``, or an unassociated
       pointer if no model with that name exists.
   * - ``print_models(self)``
     - Write the names of all registered models to standard output.

C–Fortran Bridge (``nuclear_model_interface``)
----------------------------------------------

*Defined in* ``src/Achilles/fortran/nuclear_model_interface.f90``

This module exposes ISO C-bound procedures consumed by the C++ ``FortranModel`` class.
Users typically do not call these directly, but must call ``register_all`` (via the C
entry point ``RegisterAll``) to make their models visible to Achilles.

The module maintains a global ``ModelFactory`` instance (``factory``) and a dynamic
vector of live model instances (``models``). Each live model is identified by a
1-based ``size_t`` index that is returned by ``CreateModel`` and passed to every
subsequent call.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - C entry point
     - Description
   * - ``RegisterAll()``
     - Register every built-in Fortran nuclear model with the factory and call
       ``factory%init()``. Called once at program startup by
       ``FortranModel::RegisterModels()``. To add a new model, add a
       ``call factory%register_model(...)`` line in this subroutine.
   * - ``ListModels()``
     - Print the names of all registered models; used by
       ``FortranModel::DisplayModels()``.
   * - ``CreateModel(name, idx)``
     - Create a model by name. Returns ``.true.`` and sets ``idx`` to the model's
       index on success; returns ``.false.`` if the name is not registered.
   * - ``InitModel(name, params, idx)``
     - Initialise model ``idx`` by calling its ``init`` procedure with config file
       path ``name`` and parameter map ``params``. Returns ``.true.`` on success.
   * - ``CleanUpModel(idx)``
     - Deallocate and nullify model ``idx``.
   * - ``CleanUpEvent(cur, size)``
     - Deallocate a dynamically-allocated current array of length ``size``.
   * - ``GetMode(idx)``
     - Return the nuclear interaction mode for model ``idx`` as a C ``int``.
   * - ``GetFrame(idx)``
     - Return the preferred reference frame for model ``idx`` as a C ``int``.
   * - ``ModelName(idx)``
     - Return a heap-allocated C string containing the model name for model ``idx``.
   * - ``ModelPS(idx)``
     - Return a heap-allocated C string containing the phase-space generator name.
   * - ``GetCurrents(idx, ...)``
     - Compute the hadronic currents for model ``idx``. Momenta are passed as a flat
       C row-major array and converted to ``fourvector`` objects before dispatching to
       the model's ``currents`` procedure. The output array ``cur`` has shape
       ``(nlorentz, nspin)`` in Fortran (column-major) which corresponds to
       ``(nspin, nlorentz)`` as seen by C++.
   * - ``GetInitialStateWeight(idx, ...)``
     - Return the initial-state weight from model ``idx``.
   * - ``GetInspireHEP(idx)``
     - Return a heap-allocated C string with the InspireHEP record identifier.

The ``FortranModel`` C++ Class
------------------------------

*Declared in* ``include/Achilles/fortran/FNuclearModel.hh``,
*implemented in* ``src/Achilles/fortran/FNuclearModel.cc``

``FortranModel`` is a concrete subclass of ``NuclearModel`` that delegates all
physics to a Fortran model via the C-bound procedures above. It is registered with the
C++ nuclear model factory under the name ``"FortranModel"``.

**Constructor**

The constructor reads two keys from the ``NuclearModel`` section of the YAML run card:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Name``
     - The name of the Fortran model to instantiate (must match one of the names passed
       to ``factory%register_model``).
   * - ``Ward``
     - The Ward–Takahashi identity correction to apply to the computed currents. Accepted
       values are ``None``, ``Coulomb``, ``Weyl``, and ``Landau``.
   * - ``ConfigFile``
     - Path to the model configuration file passed verbatim to the Fortran ``init``
       procedure.
   * - ``ModelParamsFile`` *(optional)*
     - Path to a YAML file whose ``Parameters`` mapping provides additional numeric
       parameters forwarded to the Fortran model.

**Overridden virtual methods**

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Behaviour
   * - ``Mode()``
     - Delegates to Fortran ``GetMode``.
   * - ``Frame()``
     - Delegates to Fortran ``GetFrame``.
   * - ``PhaseSpace(nuc_pid)``
     - Returns ``"Coherent"`` for hydrogen and free-neutron targets (setting internal
       flags used by ``CalcCurrents``); otherwise delegates to Fortran ``ModelPS``.
   * - ``CalcCurrents(...)``
     - Evaluates form factors, packs particle data into flat arrays, calls Fortran
       ``GetCurrents``, and applies the requested Ward gauge correction.
   * - ``InitialStateWeight(...)``
     - Delegates to Fortran ``GetInitialStateWeight``.
   * - ``GetName()``
     - Delegates to Fortran ``ModelName``.
   * - ``PSName()``
     - Delegates to Fortran ``ModelPS``.
   * - ``InspireHEP()``
     - Delegates to Fortran ``GetInspireHEP``.

Ward Gauge Corrections
~~~~~~~~~~~~~~~~~~~~~~

After computing the hadronic current the ``FortranModel`` can enforce a gauge condition
on the Lorentz components. The choice is controlled by the ``Ward`` run-card key:

- ``None`` — no correction applied.
- ``Coulomb`` — set the longitudinal (time-like) component to satisfy current conservation
  in the Coulomb gauge.
- ``Weyl`` — apply the Weyl (temporal) gauge condition.
- ``Landau`` — project out the longitudinal part in the Landau gauge.

Implementing a New Fortran Model
---------------------------------

To add a new nuclear model written in Fortran, follow these steps:

1. **Create a Fortran module** that ``use``\s ``nuclear_model`` and defines a concrete
   type that ``extends(model)`` and implements all deferred procedures described in
   :ref:`fortran-interface`.

   .. code-block:: fortran

      module my_model_mod
          use nuclear_model
          implicit none

          type, extends(model) :: my_model
          contains
              procedure :: init       => my_init
              procedure :: cleanup    => my_cleanup
              procedure :: mode       => my_mode
              procedure :: model_name => my_name
              procedure :: ps_name    => my_psname
              procedure :: currents   => my_currents
              procedure :: init_wgt   => my_init_wgt
              procedure :: inspirehep => my_inspirehep
          end type my_model

      contains
          ! ... implementations ...

          function build_my_model() result(ptr)
              class(model), pointer :: ptr
              allocate(my_model :: ptr)
          end function

      end module my_model_mod

2. **Register the model** in ``nuclear_model_interface.f90``. Add a ``use`` for your
   new module and a ``call factory%register_model(...)`` line inside ``register_all``:

   .. code-block:: fortran

      use my_model_mod
      ! inside register_all():
      call factory%register_model("MyModel", build_my_model)

3. **Add the source file** to the CMake build so it is compiled and linked.

4. **Reference the model** in the run card:

   .. code-block:: yaml

      NuclearModel:
        Type: FortranModel
        Name: MyModel
        Ward: None
        ConfigFile: /path/to/config.dat

.. _sherpa-interface:

Sherpa Interface
================
