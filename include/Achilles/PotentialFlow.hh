// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef POTENTIAL_FLOW_HH
#define POTENTIAL_FLOW_HH

#include <cmath>
#include <complex>

#include "Achilles/Potential.hh"
#include "Achilles/ThreeVector.hh"

namespace achilles::flow {

/// Hamiltonian of a particle of mass m in the scalar-vector potential of the nucleus,
/// H(q, p) = Re sqrt(p^2 + (m + S + i W_S)^2) + V, with S, V evaluated at (|p|, |q|).
/// This is the form that the cascade integrator differentiates; Potential::Hamiltonian of
/// the derived (non-relativistic) potentials differs from it (see the study note).
inline double Hamiltonian(const ThreeVector &q, const ThreeVector &p, const Potential &pot,
                          double mass) {
    const auto vals = pot(p.P(), q.P());
    const auto mass_eff = mass + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
    return std::sqrt(p.P2() + mass_eff * mass_eff).real() + vals.rvector;
}

/// dH/dp: the velocity. Central, momentum-dependent potential: along p-hat.
inline ThreeVector dHdp(const ThreeVector &q, const ThreeVector &p, const Potential &pot,
                        double mass) {
    const double pmag = p.P();
    if(pmag <= 0) return {0, 0, 0};
    const auto vals = pot(pmag, q.P());
    const auto dpot_dp = pot.derivative_p(pmag, q.P());
    const auto mass_eff = mass + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
    const double numerator = (vals.rscalar + mass) * dpot_dp.rscalar + pmag;
    const double denominator = std::sqrt(mass_eff * mass_eff + p.P2()).real();
    return (numerator / denominator + dpot_dp.rvector) * p / pmag;
}

/// dH/dr: minus the force. Central potential: along r-hat (not p-hat), built from
/// derivative_r. Vanishes at the origin, where a smooth central potential is flat.
inline ThreeVector dHdr(const ThreeVector &q, const ThreeVector &p, const Potential &pot,
                        double mass) {
    const double r = q.P();
    if(r <= 0) return {0, 0, 0};
    const auto vals = pot(p.P(), r);
    const auto dpot_dr = pot.derivative_r(p.P(), r);
    const auto mass_eff = mass + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
    const double numerator = (vals.rscalar + mass) * dpot_dr.rscalar;
    const double denominator = std::sqrt(mass_eff * mass_eff + p.P2()).real();
    return (numerator / denominator + dpot_dr.rvector) * q / r;
}

} // namespace achilles::flow

#endif
