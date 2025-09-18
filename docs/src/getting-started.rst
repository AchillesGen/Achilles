.. _Getting Started:

###############
Getting Started
###############


.. contents::
   :local:

.. _Obtaining the code:

******************
Obtaining the code
******************


.. _Installation:

************
Installation
************

.. _Building Achilles:

=================
Buidling Achilles
=================

In this section, all the details for building the Achilles code is provided.
The Achilles code uses CMake as a means to provide a platform agnositic installation procedure.
If you encounter any issues during the build process, please file a `github issue <https://github.com/AchillesGen/Achilles/issues>`_.

In order to build Achilles, a C++ compiler with support for `C++17 <https://en.cppreference.com/w/cpp/compiler_support/17.html>`_ is required.
In addition to a C++ compiler, a Fortran90 compiler is required (we recommend using gfortran).
Finally, Achilles requires `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ for parsing some of the data files.
All additional required dependencies will be automatically downloaded and installed for you.
If you wish to install them yourself, please see :ref:`All Dependencies` for details.

In addition to the required dependencies, Achilles has some additional optional dependencies to extend the functionality of the code.
If the user wishes to study Beyond the Standard Model processes or the decay of tau neutrinos, an additional dependency on `Sherpa`_ is required.
The interface and current limitations are described in detail in :cite:`Isaacson:2021xty` for BSM and :cite:`Isaacson:2023gwp` for tau decays.
Additional details on the interface can be found in the :ref:`Sherpa Interface <sherpa-interface>` section.
Details on enabling the interface and options for finding the Sherpa installation can be found in :ref:`Build Options`.

To build Achilles with the default options can be accomplished with the standard CMake pipeline, assuming you are in the Achilles root directory:

.. code-block:: shell-session

   $ cmake -S . -B build <additional options>
   $ cmake --build build -- -jN

Where in the above, the additional options can again be found in the :ref:`Build Options` section, and the command ``-jN`` tells CMake how many parallel jobs to be launched during the build process.
While the above commands will completely build Achilles, it will not do so in a location independent way.
To obtain this, the Achilles code needs to be installed with

.. code-block:: shell-session

   $ cmake --install build

Then Achilles will be installed into the variable set by ``CMAKE_PREFIX_PATH``, which defaults to ``/usr/local/`` on UNIX machines.
Achilles can then be launched anywhere as long as the executable can be found in your ``PATH``.


.. _All Dependencies:

============
Dependencies
============

---------------------
Required Dependencies
---------------------


---------------------
Optional Dependencies
---------------------

^^^^^^
Sherpa
^^^^^^


.. _Build Options:

=============
Build Options
=============

.. dropdown:: Achilles Options
   :open:

    +------------------------------------+--------------------------------------------+----------+
    |   Option                           |  Meaning                                   | Default  |
    +====================================+============================================+==========+
    | ``ACHILLES_ENABLE_TESTING``        | Build the Achilles test suite              | OFF      |
    +------------------------------------+--------------------------------------------+----------+
    | ``ACHILLES_ENABLE_CASCADE_TEST``   | Build the Achilles cascade executable      | OFF      |
    |                                    | to run hadron-nucleus interactions         |          |
    |                                    | or transparency checks                     |          |
    +------------------------------------+--------------------------------------------+----------+
    | ``ACHILLES_ENABLE_POTENTIAL_TEST`` | Build the Achilles potential               | OFF      |
    |                                    | executable used to test different          |          |
    |                                    | nuclear potentials in the cascade          |          |
    +------------------------------------+--------------------------------------------+----------+
    | ``ACHILLES_ENABLE_PRECOMPUTED``    | Build the Achilles interface               | OFF      |
    |                                    | to pass in pre-computed events in          |          |
    |                                    | order to cascade them using the            |          |
    |                                    | Achilles cascade                           |          |
    +------------------------------------+--------------------------------------------+----------+
    | ``ACHILLES_ENABLE_SHERPA``         | Build the Sherpa interface to              | OFF      |
    |                                    | Achilles. This requires an                 |          |
    |                                    | installation of Sherpa. Provides           |          |
    |                                    | the ability for tau decays and BSM.        |          |
    |                                    | For details see                            |          |
    |                                    | :ref:`Sherpa Interface <sherpa-interface>` |          |
    +------------------------------------+--------------------------------------------+----------+
    | ``ACHILLES_BUILD_DOCS``            | Build the Achilles manual. Requires        | OFF      |
    |                                    | doxygen, sphinx, and ...                   |          |
    +------------------------------------+--------------------------------------------+----------+
    | ``ACHILLES_EVENT_DETAILS``         | Provides additional event details useful   | OFF      |
    |                                    | for in-depth debugging.                    |          |
    +------------------------------------+--------------------------------------------+----------+

.. dropdown:: CMake Options

    +------------------------------------+--------------------------------------------+-----------------+
    |   Option                           |  Meaning                                   | Default         |
    +====================================+============================================+=================+
    | ``CMAKE_BUILD_TYPE``               | Whether to build a Release, Debug,         | Release         |
    |                                    | or RelWithDebInfo version of the code.     |                 |
    +------------------------------------+--------------------------------------------+-----------------+
    | ``CMAKE_INSTALL_PREFIX``           | Specifies the directory where the Achilles | ``/usr/local``  |
    |                                    | will be installed when ``make install``    |                 |
    |                                    | is executed                                |                 |
    +------------------------------------+--------------------------------------------+-----------------+

.. _Running:

*******
Running
*******

The main Achilles executable for lepton-nucleus interactions can be run by
launching the ``achilles`` executable, which can either be found in the
``bin`` directory of the build directory, or in the ``bin`` directory of the install directory.
Launching the code with no command line arguments, attempts to find a run card file named ``run.yml`` that
specifies the setup for the desired run. For a detailed description on the run card can be found
:ref:`here <Run Card>`.

If you wish to run with a different run card, then Achilles can be launched as

.. code-block:: shell-session
   $ ./achilles <run_card>

In which the ``<run_card>`` is the desired `YAML <https://yaml.org/>`_ configuration file to be used.

.. note:: Alternate Run Modes
   To run Achilles in other run mode configurations please see :ref:`Alternate Run Modes`.

.. dropdown:: Runtime Options
   :open:

    +--------------------------+--------------------------------+
    |   Option                 | Meaning                        |
    +==========================+================================+
    | ``-h --help``            | Show list of options           | 
    +--------------------------+--------------------------------+
    | ``--version``            | Show the version of Achilles   | 
    +--------------------------+--------------------------------+
    | ``-v`` / ``-vv``         | Increase verbosity level       | 
    +--------------------------+--------------------------------+
    | ``-l`` / ``-ll``         | Increase log verbosity level   |
    +--------------------------+--------------------------------+
    | ``--logfile``            | File to write log info to.     |
    |                          | Defaults to ``achilles.log``   |
    +--------------------------+--------------------------------+
    | ``-s --sherpa``          | Option to pass to Sherpa       |
    +--------------------------+--------------------------------+
    | ``--display-cuts``       | List available cut options     |
    +--------------------------+--------------------------------+
    | ``--display-ps``         | List available phase spaces    |
    +--------------------------+--------------------------------+
    | ``--display-ff``         | List available form factors    |
    +--------------------------+--------------------------------+
    | ``--display-int-models`` | List available cascade models  |
    |                          | for cross section evaluations  |
    +--------------------------+--------------------------------+
    | ``--display-nuc-models`` | List available nuclear models  |
    +--------------------------+--------------------------------+



.. _Dockerized Achilles:

*******************
Dockerized Achilles
*******************



.. _Sherpa: https://sherpa-team.gitlab.io/
