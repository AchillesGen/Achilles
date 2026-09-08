<!--
SPDX-FileCopyrightText: 2018-2026 Achilles Developers
SPDX-License-Identifier: CC0-1.0
-->

# Third-party components

Achilles incorporates material from the projects below.  Each file carries its
own SPDX headers; this page records the provenance in one place.

## GENIE -- pion-nucleus optical potential

`include/Achilles/OsetCrossSections.hh`, `src/Achilles/OsetCrossSections.cc`

Adapted from GENIE's `INukeOsetFormula`, written by Tomasz Golan (2015).
Copyright (c) 2003-2025, The GENIE Collaboration.  Distributed under GPLv3 with
the MCnet Guidelines for Fair Academic Usage.
<https://github.com/GENIE-MC/Generator>

The Oset parametrisation coefficients are published results of Oset et al. (see
`CITATIONS.bib`) rather than GENIE's contribution.  What follows GENIE's
implementation is the code: the arrangement of the Delta self-energy and
propagator, the delta-reduction factor, and this file's naming and commentary.
Modified by the Achilles Developers: average-nucleon kinematics
(<p^2> = 0.6 k_F^2) in place of the nucleon-at-rest approximation, the
nine-channel quasi-elastic isospin map, s-wave absorption, and integration with
the Achilles `Event` and `Nucleus` interfaces.

## ANL-Osaka -- dynamical coupled-channels amplitudes

`src/Achilles/fortran/amp_dcc_sl.f`, `src/Achilles/fortran/currents_pi_dcc.f90`

The ANL-Osaka DCC amplitude code, by S. X. Nakamura, H. Kamano and T. Sato,
provided to the Achilles Developers by its authors and distributed here with
their permission under the terms of this project's licence.

## gzstream

`src/gzstream/` -- Copyright (c) 2001 Deepak Bandyopadhyay, Lutz Kettner.
LGPL-2.1-or-later.

## Physics models implemented from the literature

Implemented from published descriptions rather than from third-party source:
the GiBUU NN elastic parametrisation; the Dmitriev-Sushkov NN -> N Delta model
(Nucl. Phys. A459 (1986) 503); the Kelly and z-expansion nucleon form factors.

See `CITATIONS.bib` for the references.
