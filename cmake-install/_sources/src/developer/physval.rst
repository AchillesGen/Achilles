.. _Physval:

##########################
Physics Validation
##########################

Physics validation (physval) checks whether a branch changes Achilles' predictions.
It generates events for a suite of experimental setups, bins them with NUISANCE3, and
compares the result with both the stored ``main`` baseline and the published data. The
summary is posted as a comment on the branch's open pull request.

This page covers day-to-day use. The statistics, the sample-specific corrections and the
storage layout are described in the `physval README`_.

.. _physval README: https://github.com/AchillesGen/Achilles/blob/main/validation/physval/README.md

.. contents::
   :local:

***************
Starting a run
***************

Put a marker in the subject line of a commit and push it to a branch. Markers in the commit
body are ignored, so a message can describe the syntax without starting a run. Markers are
also ignored on ``main``, including those a merge brings in.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Commit subject contains
     - What runs
   * - ``!physval``
     - The whole suite
   * - ``!physval(A,B)``
     - Only setups ``A`` and ``B``
   * - ``!physval(dry-run)``
     - The whole suite with synthetic events, to test the pipeline itself
   * - ``!physval(dry-run,A)``
     - Setup ``A`` with synthetic events

Setup names are the ``name`` entries in ``validation/physval/measurements.yml``. An
unknown name fails the run and lists the valid ones.

A run can also be started by hand from the *Physics Validation (NUISANCE3)* workflow on
the repository's Actions tab, which additionally lets you set the number of events and
the seed.

*******************
Reading the comment
*******************

Each measurement gets one row, and a row is flagged when main and the branch differ
significantly after correcting for the number of measurements compared.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Mark
     - Meaning
   * - 🚩
     - Changed, and agrees worse with the data
   * - ⭐
     - Changed, and agrees better with the data
   * - 🔀
     - Changed, with little effect on the agreement with the data (:math:`|\Delta\chi^2| < 1`)
   * - ✅
     - Compatible with main
   * - ⚪
     - No usable baseline, so it is compared with the data only

The columns are:

``p_cmp``
   The compatibility of main and the branch.

``q``
   ``p_cmp`` after a Benjamini-Hochberg correction over the suite. A row is flagged when
   :math:`q < 0.05`, so about 5% of flags are expected to be noise.

``Δχ²``
   The change in agreement with the data. Positive means the branch agrees worse. It only
   labels a flag and never causes one.

``p_data``
   The branch's goodness of fit to the data, for context.

A setup that fails to generate events is listed under *Did not run* with the reason, and
the workflow run fails. So does a measurement that selects no events, listed under *No
events selected*. If a sample could plausibly be empty by chance, give the reason under
``may_be_empty`` in ``measurements.yml`` and it becomes a warning.

*******************
Baselines
*******************

A nightly run on ``main`` regenerates the baselines when the physval image, the physics
code or ``measurements.yml`` has changed. A new image starts with no baseline, so
comparisons are against the data only until the next nightly run. To make one sooner, start the workflow by hand on
``main`` with the mode set to ``baseline``.

Other changes under ``validation/physval`` do not regenerate the baselines. If one moves
the predictions, start the workflow by hand in the same way, which always regenerates.

*******************
Running locally
*******************

The tests and a synthetic run need only a few Python packages:

.. code-block:: shell-session

   $ cd validation/physval
   $ pip install numpy scipy pyyaml matplotlib pytest
   $ python3 -m pytest
   $ python3 physval.py --dry-run --make-baseline --key dev --out-dir /tmp/pv/baselines/dev
   $ python3 physval.py --dry-run --key dev --baseline-dir /tmp/pv --out-dir /tmp/pv/out

A real run needs ``achilles`` on the ``PATH`` and pyNUISANCE, which are available in the
``ghcr.io/achillesgen/achilles-physval`` image.
