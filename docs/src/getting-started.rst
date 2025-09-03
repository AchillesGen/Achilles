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


.. _All Dependencies:

----------------
All Dependencies
----------------


.. _Build Options:

-------------
Build Options
-------------


.. _Dockerized Achilles:

===================
Dockerized Achilles
===================


.. _Running:

*******
Running
*******

.. _Sherpa: https://sherpa-team.gitlab.io/
