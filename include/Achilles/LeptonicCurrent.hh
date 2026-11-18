// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "Achilles/Current.hh"

#include <complex>
#include <map>
#include <vector>

namespace achilles {

template <int D> class FourVectorT;
using FourVector = FourVectorT<1>;
class Particle;
class Nucleus;
class Event;
class NuclearModel;
class PID;
struct FormFactorInfo;
struct ProcessInfo;

using Particles = std::vector<Particle>;
using Current = std::vector<VCurrent>;
using Currents = std::map<int, Current>;
using FFDictionary = std::map<std::pair<PID, PID>, std::vector<FormFactorInfo>>;

class LeptonicCurrent {
  public:
    LeptonicCurrent() = default;
    void Initialize(const ProcessInfo &);
    FFDictionary GetFormFactor();
    Currents CalcCurrents(const FourVector &, const FourVector &) const;
    size_t NSpins() const { return 4; }

  private:
    bool NeutralCurrent(PID, PID) const;
    bool ChargedCurrent(bool, PID, PID) const;
    std::complex<double> coupl_left{}, coupl_right{};
    double mass{}, width{};
    int pid{};
    bool anti{};
};

} // namespace achilles
