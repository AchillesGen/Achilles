#include "Achilles/MomSolver.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Potential.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

#include "spdlog/spdlog.h"

achilles::FourVector achilles::SolveDelta(const achilles::FourVector &p1,
                                          const achilles::FourVector &p2, double m3, double m4,
                                          double cosTheta, double phi) {
    const double s = (p1 + p2).M2();
    const auto beta = (p1 + p2).BoostVector();
    const auto pcm = (p1 + p2).Boost(-beta);
    const double sinTheta = sqrt(1 - pow(cosTheta, 2));
    auto _func = [&](double energy) {
        double p3Mag = sqrt(energy * energy - m3 * m3);
        FourVector p3(energy, p3Mag * sinTheta * cos(phi), p3Mag * sinTheta * sin(phi),
                      p3Mag * cosTheta);
        return s - 2 * pcm * p3 + m3 * m3 - m4 * m4;
    };

    Brent brent(_func);
    double energy = brent.CalcRoot(m3, (p1 + p2).E());
    double p3Mag = sqrt(energy * energy - m3 * m3);
    return FourVector(energy, p3Mag * sinTheta * cos(phi), p3Mag * sinTheta * sin(phi),
                      p3Mag * cosTheta)
        .Boost(beta);
}

std::pair<double, double> achilles::FindMomentumRange(const achilles::FourVector &q_free,
                                                      double rangeExtend) {
    std::pair<double, double> range;
    const double detM = sqrt((q_free.M2() - 4 * pow(achilles::Constant::mN, 2)) *
                             (q_free.M2() + q_free.P2()) / q_free.M2());
    range.first = q_free.M2() > 2 * achilles::Constant::mN * q_free.E() ? (detM - q_free.P()) / 2
                                                                        : (q_free.P() - detM) / 2;
    range.second = (q_free.P() + detM) / 2;
    range.first /= rangeExtend;
    range.second *= rangeExtend;

    return range;
}

achilles::FourVector achilles::SolveDeltaWithPotential(const achilles::FourVector &q,
                                                       const Potential &potential, double m3,
                                                       double m4, double p3Mag, double phi,
                                                       double radius3, double radius4) {
    auto potential3 = potential(p3Mag, radius3);
    const double energy =
        sqrt(p3Mag * p3Mag + pow(m3 + potential3.rscalar, 2)) + potential3.rvector;
    auto _func = [&](double cosTheta) {
        double sinTheta = sqrt(1 - cosTheta * cosTheta);
        FourVector p3 = {energy, p3Mag * sinTheta * cos(phi), p3Mag * sinTheta * sin(phi),
                         p3Mag * cosTheta};
        auto p4 = q - p3;
        auto potential4 = potential(p4.P(), radius4);
        p4.E() = sqrt(p4.P2() + pow(m4 + potential4.rscalar, 2)) + potential4.rvector;
        return q.E() - p3.E() - p4.E();
    };

    Brent brent(_func);
    double cosTheta = brent.CalcRoot(-1, 1);
    double sinTheta = sqrt(1 - cosTheta * cosTheta);
    return {energy, p3Mag * sinTheta * cos(phi), p3Mag * sinTheta * sin(phi), p3Mag * cosTheta};
}
