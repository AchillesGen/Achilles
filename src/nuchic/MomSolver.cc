#include "nuchic/MomSolver.hh"
#include "nuchic/Constants.hh"
#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/Utilities.hh"

#include "spdlog/spdlog.h"

nuchic::FourVector nuchic::SolveDelta(const nuchic::FourVector &p1, const nuchic::FourVector &p2,
                                      double m3, double m4, double cosTheta, double phi) {
    const double s = (p1 + p2).M2(); 
    const auto beta = (p1 + p2).BoostVector();
    const auto pcm = (p1 + p2).Boost(-beta);
    const double sinTheta = sqrt(1 - pow(cosTheta, 2));
    auto _func = [&](double energy) {
        double p3Mag = sqrt(energy*energy - m3*m3);
        FourVector p3(p3Mag*sinTheta*cos(phi),
                p3Mag*sinTheta*sin(phi), p3Mag*cosTheta, energy);
        return s - 2*pcm*p3 + m3*m3 - m4*m4;
    };

    Brent brent(_func);
    double energy = brent.CalcRoot(m3, (p1 + p2).E());
    double p3Mag = sqrt(energy*energy - m3*m3);
    return FourVector(p3Mag*sinTheta*cos(phi), p3Mag*sinTheta*sin(phi), p3Mag*cosTheta, energy).Boost(beta);
}

double nuchic::Potential(double p, double rho) {
    constexpr double rho0 = 0.16;
    rho = rho/rho0;
    const double alpha = 15.52 * rho + 24.93*rho*rho;
    const double beta = -116*rho;
    const double lambda = 3.29 - 0.373*rho;
    const double k = p/nuchic::Constant::HBARC;
    return alpha + beta/(1+pow(k/lambda, 2));
}

std::pair<double, double> nuchic::FindMomentumRange(double E, double q, double m,
                                                    double rho1, double rho2) {
    auto _func = [&](double p3, double sign) {
        return E - sqrt(m*m+p3*p3) - sqrt(m*m+q*q-2*sign*q*p3+p3*p3)
            - Potential(p3, rho1) - Potential(sqrt(q*q-2*sign*q*p3+p3*p3), rho2);
    };

    auto _funcMin = [&](double p3) {
        return _func(p3, -1);
    };

    auto _funcMid = [&](double p3) {
        return _func(p3, 0);
    };

    auto _funcMax = [&](double p3) {
        return _func(p3, 1);
    };

    spdlog::info("Min: {}, {}", _funcMin(0), _funcMin(2*q));
    spdlog::info("Mid: {}, {}", _funcMid(0), _funcMid(2*q));
    spdlog::info("Max: {}, {}", _funcMax(0), _funcMax(2*q));
    for(size_t i = 0; i < 100; ++i) {
        double tmp = static_cast<double>(i)*2*q/100;
        spdlog::info("_funcMin({}) = {}", tmp, _funcMin(tmp));
    }
    for(size_t i = 0; i < 100; ++i) {
        double tmp = static_cast<double>(i)*2*q/100;
        spdlog::info("_funcMid({}) = {}", tmp, _funcMid(tmp));
    }
    for(size_t i = 0; i < 100; ++i) {
        double tmp = static_cast<double>(i)*2*q/100;
        spdlog::info("_funcMax({}) = {}", tmp, _funcMax(tmp));
    }
    Brent brentMin(_funcMin);
    double pmin = brentMin.CalcRoot(0, 2*q);
    spdlog::info("Min sol: {}", pmin);
    Brent brentMax(_funcMax);
    double pmax = brentMax.CalcRoot(0, 2*q);
    
    return {pmin, pmax};
}

nuchic::FourVector nuchic::SolveDeltaWithPotential(const nuchic::FourVector &q,
                                                   double m3, double m4,
                                                   double p3Mag, double phi,
                                                   double rho3, double rho4) {
    const double energy = sqrt(p3Mag*p3Mag + m3*m3)+Potential(p3Mag, rho3);
    auto _func = [&](double cosTheta) {
        double sinTheta = sqrt(1-cosTheta*cosTheta);
        FourVector p3 = {p3Mag*sinTheta*cos(phi), p3Mag*sinTheta*sin(phi), p3Mag*cosTheta, energy};
        auto p4 = q - p3;
        p4.E() = sqrt(p4.P2() + m4*m4) + Potential(p4.P(), rho4);
        return q.E() - p3.E() - p4.E();
    };

    Brent brent(_func);
    double cosTheta = brent.CalcRoot(-1, 1);
    double sinTheta = sqrt(1-cosTheta*cosTheta);
    return {p3Mag*sinTheta*cos(phi), p3Mag*sinTheta*sin(phi), p3Mag*cosTheta, energy};
}
