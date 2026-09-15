<!--
SPDX-FileCopyrightText: 2018-2026 Achilles Developers
SPDX-License-Identifier: CC0-1.0
-->

# Third-party components

Achilles incorporates material from the projects below.  Each file carries its
own SPDX headers; this page records the provenance in one place.

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
(Nucl. Phys. A459 (1986) 503); the Oset et al. pion absorption and quasielastic
model (Nucl. Phys. A468 (1987) 631; A484 (1988) 557); the Kelly and z-expansion
nucleon form factors.

See `CITATIONS.bib` for the references.
