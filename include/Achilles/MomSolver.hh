// SPDX-FileCopyrightText: 2018-2024 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MOM_SOLVER_HH
#define MOM_SOLVER_HH

#include <memory>
#include <utility>

namespace achilles {
class FourVector;
class Potential;
class Nucleus;

FourVector SolveDelta(const FourVector &, const FourVector &, double, double, double, double);
std::pair<double, double> FindMomentumRange(const FourVector &, double rangeExtend = 1.05);
FourVector SolveDeltaWithPotential(const FourVector &, const Potential &, double, double, double,
                                   double, double, double);
} // namespace achilles

#endif
