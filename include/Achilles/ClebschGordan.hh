// SPDX-FileCopyrightText: 2018-2024 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef ACHILLES_CLEBSCHGORDAN
#define ACHILLES_CLEBSCHGORDAN

namespace achilles {

struct SpinState {
    int total, zaxis;
    SpinState() = default;
    SpinState(double _total, double _zaxis)
        : total{static_cast<int>(2 * _total)}, zaxis{static_cast<int>(2 * _zaxis)} {}
};

double ClebschGordan(SpinState, SpinState, SpinState);

} // namespace achilles

#endif
