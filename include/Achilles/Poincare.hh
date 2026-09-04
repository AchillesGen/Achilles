// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "Achilles/FourVector.hh"

#include <array>
#include <vector>

namespace achilles {

/// Poincare transformations are pure geometry: depending on the mode its
/// reference vectors are either momenta or dimensionless directions, so it
/// works on canonical component values (MeV) and converts at its FourVector
/// interface rather than carrying a mass dimension of its own.
class Poincare {
  private:
    using Vec4 = std::array<double, 4>;

    int m_type;
    Vec4 m_l{}, m_t{};
    double m_rsq, m_omct{}, m_st{};

  public:
    Poincare(const FourVector &v = FourVector(), const double &rsq = -1.);
    Poincare(const FourVector &v1, const FourVector &v2, int mode = 0);

    void Boost(FourVector &v) const;
    void BoostBack(FourVector &v) const;

    void Rotate(FourVector &v) const;
    void RotateBack(FourVector &v) const;

    void Lambda(FourVector &v) const;
    void LambdaBack(FourVector &v) const;

    void Invert();

    inline FourVector operator*(const FourVector &vin) const {
        FourVector v(vin);
        if(m_type == 1) Boost(v);
        if(m_type == 2) Rotate(v);
        if(m_type == 3) Lambda(v);
        return v;
    }

    inline FourVector PL() const { return FourVector::FromNative(m_l); }
    inline FourVector PT() const { return FourVector::FromNative(m_t); }

    inline double OMCTheta() const { return m_omct; }
    inline double SinTheta() const { return m_st; }

}; // end of class Poincare

class Poincare_Sequence : public std::vector<Poincare> {
  public:
    FourVector operator*(const FourVector &p) const;
    void Invert();

}; // end of class Poincare_Sequence

} // namespace achilles
