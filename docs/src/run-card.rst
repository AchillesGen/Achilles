.. _Run Card Structure:

##################
Run Card Structure
##################

An Achilles simulation is steered by a set of options. These options enable the user to run Achilles
in various ways to study the different physics choices in the code.

The options have to be specified in a run card which by default is named ``run.yml`` residing in the current
working directory. If you want to use a different run card, it can be passed into the Achilles executable
as the first positional argument (``achilles <input_file>``).

All of the configurations options are written in the `YAML <https://yaml.org/>`_ format. The different
sections in the run card are described below, along with the currently available options, and references
to the physics where appropriate.

.. warning::
   Achilles is still in beta development, and the options are subject to change before the v1.0 release
   of Achilles.

.. contents::
   :local:


***********************
Including External Files
***********************

Any value in the run card can be replaced by the contents of a separate YAML file using the
custom ``!include`` tag:

.. code-block:: yaml

   SomeKey: !include "path/to/file.yml"

When Achilles loads the run card it replaces every ``!include`` scalar with the fully parsed
contents of the referenced file. The path is resolved relative to the Achilles data search
paths (first in the current working directory, then the environment variable ``ACHILLES_PATH``, and finally the install location); absolute paths are also accepted. Included files are
processed recursively, so an included file may itself contain further ``!include`` tags.

The ``!include`` tag is particularly useful for sharing nucleus definitions, cascade
interaction tables, form factor files, and option defaults across multiple run cards without
duplicating their contents. The pre-supplied defaults under ``data/default/`` make
extensive use of this mechanism:

.. code-block:: yaml

   Nuclei:
     - Nucleus: !include "data/default/12C.yml"

   Cascade:
     Interactions: !include "data/default/VirtResInteractions.yml"

   Options: !include "data/default/OptionDefaults.yml"


*************
Main Settings
*************

The first section in the run card is the ``Main`` section. This section contains information about the basic options for running Achilles. These include:

* ``NEvents`` (integer): Number of unweighted events to generate
* ``HardCuts`` (boolean): Flag to turn on / off cuts at the hard scattering level
* ``EventCuts`` (boolean): :yellow:`[deprecated]` Handles cuts to be applied after all simulation steps
* ``RunDecays`` (boolean): Decay particles that are unstable longer than the typical cascade time scale,
  such as tau leptons (requires Sherpa to be linked)
* ``Output`` (YAML Node): Parameters determining how events are written out to file
    * ``Format`` (string): Format of file to be written out. Current supported options are:
      ``NuHepMC`` (default), ``HepMC`` (:yellow:`[deprecated]`), and ``Achilles``.
    * ``Name`` (string): File path from the current working directory to write the events out to
    * ``Zipped`` (boolean): Flag to determine if the output should be written compressed using the
      gzip library or not.


Additional discussion on the output formats are discussed in detail in the `Output Formats`_ section.

Example:

.. code-block:: yaml

   Main:
     NEvents: 50000
     HardCuts: true
     RunDecays: false
     Output:
       Format: NuHepMC
       Name: achilles.hepmc
       Zipped: true


*********
Processes
*********

The processes section is where the set of incoming and outgoing leptons are listed.
Each line in this section denotes a process to calculate. Achilles will combine this information
with the list of nuclear models to consider listed in `NuclearModels`_.

An example setup for running both neutral current and charged current is given below:

.. code-block:: yaml

   Processes:
     - Leptons: [14, [14]]
     - Leptons: [14, [13]]

The above code block can be broken down as follows. First is the YAML node giving the name of the section
in this case ``Processes``. It is followed by a YAML list item. Each list item is described first by
denoting that it is describing the leptons in the process, and then followed by a list of two elements.
The first item is the PID for the incoming lepton. The second item is a list of all outgoing leptons.

.. note::

    In order to handle processes with final states beyond single leptons, Achilles needs to be compiled
    with Sherpa support to handle the n-body final state phase space.

The allowed values for the leptons are any that are defined in the model using the standard
Particle Data Group Particle Identification Codes (`PID`_) numbering scheme.
This section is designed to be flexible enough to define arbitrary processes that the user
would be interested in studying through the use of the Achilles event generator.
The extension to Beyond the Standard Model particles is also supported in this section
through the use of the :ref:`Sherpa interface <sherpa-interface>`. For details on
defining the physics model to be used, please see the discussion in `Sherpa Options`_.

*****
Beams
*****

The ``Beams`` section defines the energy distribution(s) of the incoming beam particle(s).
Each entry in the list specifies a single beam via its PDG ``PID`` and a ``Beam Params``
sub-node that selects the beam type and its parameters.

.. code-block:: yaml

   Beams:
     - Beam:
         PID: 14
         Beam Params:
           Type: Monochromatic
           Energy: 1000

Multiple entries are allowed when running with more than one beam species simultaneously;
however, all beams must share the same energy range.

``Beam Params`` types
---------------------

**Monochromatic**

A fixed-energy beam. All events use exactly the specified energy.

.. code-block:: yaml

   Beam Params:
     Type: Monochromatic
     Energy: 1108          # beam energy in MeV

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Key
     - Description
   * - ``Energy``
     - Beam energy in MeV.

**Spectrum**

A beam whose energy is sampled from a flux histogram read from a data file. The file
format is auto-detected from its header line; supported formats are ``Achilles``,
``MiniBooNE``, and ``T2K`` (ND280).

.. code-block:: yaml

   Beam Params:
     Type: Spectrum
     Histogram: flux/my_flux.txt

Alternatively, a flux histogram stored in a HepData YAML file can be used:

.. code-block:: yaml

   Beam Params:
     Type: Spectrum
     HepData: flux/my_flux.yaml

If Achilles is compiled with ROOT support, a ROOT histogram can be specified:

.. code-block:: yaml

   Beam Params:
     Type: Spectrum
     ROOTHist:
       File: my_flux.root      # resolved under flux/ directory

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Key
     - Description
   * - ``Histogram``
     - Path to a text-format flux histogram. The file format (``Achilles``,
       ``MiniBooNE``, or ``T2K``) is determined from the file's header line.
   * - ``HepData``
     - Path to a HepData YAML file containing the flux histogram.
   * - ``ROOTHist``
     - Sub-node with a ``File`` key giving the ROOT file name (requires ROOT support).

**FlatFlux**

A beam with a uniform (flat) energy distribution between ``MinEnergy`` and ``MaxEnergy``.

.. code-block:: yaml

   Beam Params:
     Type: FlatFlux
     MinEnergy: 200          # MeV
     MaxEnergy: 2000         # MeV

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Key
     - Description
   * - ``MinEnergy``
     - Minimum beam energy in MeV.
   * - ``MaxEnergy``
     - Maximum beam energy in MeV.

*******
Cascade
*******

The ``Cascade`` section controls the intranuclear cascade (INC) that propagates the
hadronic final-state particles through the nucleus after the hard interaction.

.. code-block:: yaml

   Cascade:
     Run: false
     Interactions: !include "data/default/VirtResInteractions.yml"
     Step: 0.04
     Probability: Cylinder
     InMedium: None
     PotentialProp: false
     Algorithm: Base

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Run``
     - Boolean. Set to ``true`` to enable the intranuclear cascade; ``false`` skips
       the cascade entirely.
   * - ``Interactions``
     - List of interaction model definitions used during the cascade. Typically
       provided via ``!include`` (see below). Each entry has a ``Name`` and an
       ``Options`` sub-node.
   * - ``Step``
     - Step size for particle propagation through the nucleus in fm. Smaller values
       increase accuracy at the cost of computation time. Default: ``0.04``.
   * - ``Probability``
     - Model used to compute the interaction probability at each step. Accepted
       values: ``Cylinder``, ``Gaussian``.
   * - ``InMedium``
     - In-medium correction applied to particle kinematics during propagation.
       Accepted values: ``None``, ``NonRelativistic``.
   * - ``PotentialProp``
     - Boolean. When ``true``, the nuclear potential is used to modify particle
       trajectories during propagation.
   * - ``Algorithm``
     - Cascade algorithm. Accepted values: ``Base`` (standard step-based propagation),
       ``MFP`` (mean-free-path sampling).

Interaction Models
------------------

The ``Interactions`` list specifies which hadronic interaction models are active
during the cascade. Each entry requires a ``Name`` and an ``Options`` sub-node.

**NucleonNucleon**

Nucleon–nucleon rescattering.

.. code-block:: yaml

   - Name: NucleonNucleon
     Options:
       Mode: GiBUU
       ResonanceMode: Decay
       ExpSup: 0.0

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Description
   * - ``Mode``
     - Cross-section parametrisation. Currently only ``GiBUU`` is supported.
   * - ``ResonanceMode``
     - Treatment of resonances produced during NN scattering. ``Decay`` causes
       resonances to decay immediately; ``Propagate`` propagates them through the
       cascade and decays probabilistically according to the particle lifetime.
   * - ``ExpSup``
     - Exponential suppression factor. Set to ``0.0`` to disable.

**DeltaInteraction**

Delta resonance interactions (used together with ``ResonanceMode: Propagate``).

.. code-block:: yaml

   - Name: DeltaInteraction
     Options:
       Mode: GiBUU
       SWaveAbsorption: true
       ExpSup: 0.0

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Description
   * - ``Mode``
     - Cross-section parametrisation. Currently only ``GiBUU`` is supported.
   * - ``SWaveAbsorption``
     - Boolean. Enable s-wave pion absorption via the Delta using the Oset model.
   * - ``ExpSup``
     - Exponential suppression factor. Set to ``0.0`` to disable.

**PionInteraction**

Pion rescattering and absorption.

.. code-block:: yaml

   - Name: PionInteraction
     Options:
       HardScatter:
         Name: MesonBaryonInteraction
         Options:
       Absorption:
         Name: PionAbsorptionOneStep
         Options:

**ConstantInteraction**

Assigns a fixed cross section to specified two-body initial states. For example, can be used to give
zero-cross-section placeholders to hyperon and photon channels:

.. code-block:: yaml

   - Name: ConstantInteraction
     Options:
       InitialStates:
         - Incoming: [22, 2112]
           Outgoing:
             - Outgoing: [22, 2112]
               CrossSection: 0

Each entry under ``InitialStates`` maps a pair of incoming PIDs to one or more
outgoing channels, each with a constant ``CrossSection`` in mb.

*************
NuclearModels
*************

The ``NuclearModels`` section is a list of nuclear model configurations. Achilles
evaluates every listed model for each process defined in `Processes`_, so multiple
models can be run simultaneously within a single job.

Each entry begins with a ``NuclearModel:`` key whose sub-node contains at minimum a
``Model`` field identifying the model type, along with model-specific options described
below. The full list of available models loaded into Achilles can be displayed by running
with the option ``--display-nuc-models``.

.. code-block:: yaml

   NuclearModels:
     - NuclearModel:
         Model: QESpectral
         SpectralP: data/Spectral_Functions/pke12p_tot.data
         SpectralN: data/Spectral_Functions/pke12n_tot.data
         FormFactorFile: "FormFactors.yml"
         Ward: None
     - NuclearModel:
         Model: FortranModel
         Name: RES_Spectral_Func
         ConfigFile: data/info_C12_pke.data
         FormFactorFile: "FormFactors.yml"
         Ward: None

Common fields
-------------

The following fields are shared by all model types:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Model``
     - Model type identifier. Built-in values: ``QESpectral``, ``HyperonSpectral``,
       ``FortranModel``.
   * - ``FormFactorFile``
     - Path to a YAML file defining the form factor parametrisations. See
       ``data/default/FormFactors.yml`` for the full set of supported parametrisations
       and their parameters.
   * - ``Ward``
     - Ward–Takahashi identity gauge correction applied to the computed hadronic
       current. Accepted values: ``None``, ``Coulomb``, ``Weyl``, ``Landau``.

``QESpectral``
--------------

Quasi-elastic scattering computed with a nucleon spectral function.

.. code-block:: yaml

   NuclearModel:
     Model: QESpectral
     SpectralP: data/Spectral_Functions/pke12p_tot.data
     SpectralN: data/Spectral_Functions/pke12n_tot.data
     FormFactorFile: "FormFactors.yml"
     Ward: None

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``SpectralP``
     - Path to the proton spectral function data file.
   * - ``SpectralN``
     - Path to the neutron spectral function data file.

``HyperonSpectral``
-------------------

Hyperon production computed with a nucleon spectral function. The configuration keys
are identical to ``QESpectral``:

.. code-block:: yaml

   NuclearModel:
     Model: HyperonSpectral
     SpectralP: data/pke12_tot.data
     SpectralN: data/pke12_tot.data
     FormFactorFile: "FormFactors.yml"
     Ward: None

``FortranModel``
----------------

A nuclear model implemented in Fortran and loaded through the
:ref:`Fortran interface <fortran-interface>`.

.. code-block:: yaml

   NuclearModel:
     Model: FortranModel
     Name: RES_Spectral_Func
     ConfigFile: data/info_C12_pke.data
     FormFactorFile: "FormFactors.yml"
     Ward: None

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Name``
     - Name of the registered Fortran model to instantiate (must match the name
       passed to ``factory%register_model`` in ``nuclear_model_interface.f90``).
   * - ``ConfigFile``
     - Path to the model configuration file passed verbatim to the Fortran ``init``
       procedure.
   * - ``ModelParamsFile`` *(optional)*
     - Path to a YAML file whose ``Parameters`` mapping supplies additional numeric
       parameters to the Fortran model.

``FormFactorFile`` format
-------------------------

The form factor YAML file (e.g. ``data/default/FormFactors.yml``) has two parts.
The first part maps each form factor category to its parametrisation name:

.. code-block:: yaml

   vector: Kelly
   axial: AxialZExpansion
   coherent: Helm
   resonancevector: ResonanceVectorDummy
   resonanceaxial: ResonanceAxialDummy
   mecvector: MesonExchangeVector
   mecaxial: MesonExchangeAxial
   hyperon: Hyperon

The second part provides the parameter values for each named parametrisation. The
available parametrisations for the vector form factors include ``Kelly``, ``BBBA``,
``ArringtonHill``, and ``VectorDipole``; for the axial form factor ``AxialZExpansion``
and ``AxialDipole`` are supported; for coherent scattering ``Helm`` and ``Lovato``
are available. Each parametrisation section must supply the numeric parameters listed
in the corresponding form factor implementation.

******
Nuclei
******

The ``Nuclei`` section is a list of target nuclei to simulate. Each entry is a
``Nucleus:`` node that can be written inline or loaded from a file with ``!include``.
The pre-supplied nucleus files are located in ``data/default/``.

.. code-block:: yaml

   Nuclei:
     - Nucleus: !include "data/default/12C.yml"

Currently provided nucleus files: ``1H.yml``, ``1N.yml``, ``12C.yml``, ``16O.yml``,
``40Ar.yml``.

Nucleus fields
--------------

.. code-block:: yaml

   Name: 12C
   Binding: 8.6            # MeV
   Fermi Momentum: 225     # MeV

   Density:
     ProtonFile: data/densities/c12.prova.txt
     NeutronFile: data/densities/c12.prova.txt
     FilePotential: data/realOP_12C_EDAI.dat
     Function: configuration
     Configs: data/configurations/QMC_configs.out.gz

   FermiGas:
     Type: Local
     Correlated: false
     SRCfraction: 0.2
     LambdaSRC: 2.75
     Params: []

   Potential:
     Name: Schroedinger
     r0: 0.16
     Mode: 3

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Name``
     - Human-readable nucleus identifier (e.g. ``12C``, ``40Ar``).
   * - ``Binding``
     - Average binding energy per nucleon in MeV.
   * - ``Fermi Momentum``
     - Fermi momentum in MeV.

**Density sub-node**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``ProtonFile``
     - Path to the proton density profile data file.
   * - ``NeutronFile``
     - Path to the neutron density profile data file.
   * - ``FilePotential``
     - Path to the real optical potential data file.
   * - ``Function``
     - Density function type. Use ``configuration`` to sample from a set of
       pre-generated nucleon configurations.
   * - ``Configs``
     - Path to the nucleon configuration file (gzip-compressed). Required when
       ``Function`` is ``configuration``.

**FermiGas sub-node**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Type``
     - Fermi gas model type. Currently ``Local`` (local Fermi gas) is supported.
   * - ``Correlated``
     - Boolean. Enable short-range correlation (SRC) corrections.
   * - ``SRCfraction``
     - Fraction of nucleons assigned to the SRC component.
   * - ``LambdaSRC``
     - SRC momentum cut-off scale in GeV.
   * - ``Params``
     - List of additional model-specific parameters (may be empty).

**Potential sub-node**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Name``
     - Potential model. Accepted values: ``Schroedinger`` (Schrödinger-equivalent
       optical potential), ``Cooper`` (relativistic optical potential), ``Wiringa``
       (non-relativistic Wiringa potential).
   * - ``r0``
     - Nuclear saturation density in fm\ :sup:`-3`. Typical value: ``0.16``.
   * - ``Mode``
     - Integer mode flag passed to the potential implementation.

********
HardCuts
********

The ``HardCuts`` section defines kinematic selection cuts applied to the hard-scattering
final state before event weighting. Cuts are only active when ``HardCuts: true`` is set
in `Main Settings`_. Each entry in the list specifies one cut:

.. code-block:: yaml

   HardCuts:
     - Type: AngleTheta
       PIDs: 11
       range: [36.5, 38.5]

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Type``
     - Cut type. ``AngleTheta`` selects events where the polar angle :math:`\theta`
       (in degrees) of the specified particle falls within the given range.
   * - ``PIDs``
     - PDG PID of the particle to which the cut is applied. A single integer or a
       list of integers is accepted.
   * - ``range``
     - Two-element list ``[min, max]`` giving the allowed range for the cut variable.

For a list of all implemented cuts, a user can run Achilles with the option ``--display-cuts``.

*******
Options
*******

The ``Options`` section controls the numerical integration and unweighting settings.
It is most conveniently loaded from a shared defaults file:

.. code-block:: yaml

   Options: !include "data/default/OptionDefaults.yml"

The full inline form is:

.. code-block:: yaml

   Options:
     Initialize:
       Seed: 12345678
       Accuracy: 1.0e-2
     Unweighting:
       Name: Percentile
       percentile: 99

**Initialize sub-node**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Seed``
     - Integer seed for the random number generator. Use different seeds to obtain
       statistically independent runs. Setting to -1 results in using a different seed each run.
   * - ``Accuracy``
     - Target relative integration accuracy for the adaptive Monte Carlo integration.
       Smaller values require more function evaluations.

**Unweighting sub-node**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Name``
     - Unweighting algorithm. ``Percentile`` (default) determines the maximum weight from a
       high percentile of the weight distribution, reducing the impact of rare
       large-weight events. In this mode, Achilles returns partially unwighted events.
       Another option is ``NoUnweigheter`` set using option ``None``, which
       make Achilles return fully weighted events.
   * - ``percentile``
     - Percentile (0–100) used to set the maximum weight when ``Name`` is
       ``Percentile``. A value of ``99`` discards the top 1% of weights.

*******
Backend
*******

The ``Backend`` section selects the cross-section calculation backend.

.. code-block:: yaml

   Backend:
     Name: Default
     Options: []

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Name``
     - Backend identifier. ``Default`` uses the built-in matrix-element evaluator.
       ``Sherpa`` and ``BSM`` delegate to the Sherpa library (requires Achilles to be
       compiled with ``ACHILLES_SHERPA_INTERFACE``).
   * - ``Options``
     - Backend-specific options passed as a YAML sequence. The ``Default`` backend
       accepts an empty list.

**************
Sherpa Options
**************

When the ``Sherpa`` or ``BSM`` backend is selected, Achilles requires a top-level
``SherpaOptions`` section. This is most conveniently provided via ``!include``:

.. code-block:: yaml

   SherpaOptions: !include "data/default/SherpaOptions.yml"

The inline form is:

.. code-block:: yaml

   SherpaOptions:
     Model: SMnu
     ParamCard: parameters.dat

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Description
   * - ``Model``
     - Sherpa model identifier. ``SMnu`` is the Standard Model setup with nucleons to not mess
       with the default Standard Model setup used in collider physics within Sherpa.
       Custom UFO models are also supported; see the
       :ref:`Sherpa interface <sherpa-interface>` documentation.
   * - ``ParamCard``
     - Path to the Sherpa parameter card file that sets coupling constants and particle
       masses for the selected model.

**************
Output Formats
**************

The output format is set by the ``Format`` key in the ``Output`` sub-node of
`Main Settings`_.

``NuHepMC``
   The recommended output format. Events are written in the `NuHepMC
   <https://github.com/NuHepMC/Spec>`_ convention on top of the HepMC3 library.
   NuHepMC files are readable by standard HepMC3 analysis tools and carry
   neutrino-experiment-specific metadata. Compressed output (``Zipped: true``) uses
   gzip.

``HepMC`` *(deprecated)*
   Plain HepMC3 format without NuHepMC metadata. Retained for backward compatibility.

``Achilles``
   Native Achilles binary format. This is a compact representation intended for
   internal debugging and post-processing within the Achilles framework.


.. _PID: https://pdg.lbl.gov/2024/reviews/rpp2024-rev-monte-carlo-numbering.pdf
