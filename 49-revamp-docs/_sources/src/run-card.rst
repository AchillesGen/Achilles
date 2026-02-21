.. _Run Card Structure:

##################
Run Card Structure
##################

An Achilles setup is steered by a various set of options. These options enable the user to run Achilles
in a various ways to study the different physics choices in the code.

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


*************
Main Settings
*************

The first section in the run card is the ``Main`` section. This section contains information about the basic options for running Achilles. These include:

* ``NEvents`` (integer): Number of unweighted events to generate
* ``HardCuts`` (boolean): Flag to turn on / off cuts at the hard scattering level
* ``EventCuts`` (boolean): :yellow:`[deprecated]` Handles cuts to be applied after all simulation steps
* ``RunDecays`` (boolean): Decay particles that are unstable longer than the typical cascade time scale,
  such as tau leptons
* ``Output`` (YAML Node): Parameters determining how events are written out to file
    * ``Format`` (string): Format of file to be written out. Current supported options are:
      ``NuHepMC``, ``HepMC`` (:yellow:`[deprecated]`), and ``Achilles``.
    * ``Name`` (string): File path from the current working directory to write the events out to
    * ``Zipped`` (boolean): Flag to determine if the output should be written compressed using the
      gzip library or not.


Additional discussion on the output formats are discussed in detail in the `Output Formats`_ section.


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

*******
Cascade
*******

*************
NuclearModels
*************


******
Nuclei
******


********
HardCuts
********


*******
Options
*******


*******
Backend
*******


**************
Sherpa Options
**************


.. _PID: https://pdg.lbl.gov/2024/reviews/rpp2024-rev-monte-carlo-numbering.pdf 
