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

double Potential(double p, double rho) {
    constexpr double rho0 = 0.16;
    rho = rho/rho0;
    const double alpha = 15.52 * rho + 24.93*rho*rho;
    const double beta = -116*rho;
    const double lambda = 3.29 - 0.373*rho;
    const double k = p/nuchic::Constant::HBARC;
    return alpha + beta/(1+pow(k/lambda, 2));
}

nuchic::FourVector nuchic::SolveDeltaWithPotential(const nuchic::FourVector &p1, const nuchic::FourVector &p2,
                                                   double m3, double m4,
                                                   double cosTheta, double phi, double /*rho*/) {
    const double s = (p1 + p2).M2(); 
    // const auto beta = (p1 + p2).BoostVector();
    // const auto pcm = (p1 + p2).Boost(-beta);
    const double sinTheta = sqrt(1 - pow(cosTheta, 2));
    auto _func = [&](double p3Mag) {
        double pot3 = 0; // Potential(p3Mag, rho);
        double energy = sqrt(p3Mag*p3Mag + m3*m3);
        FourVector p3(p3Mag*sinTheta*cos(phi),
                p3Mag*sinTheta*sin(phi), p3Mag*cosTheta, energy + pot3);
        // FourVector p4 = p1 + p2 - p3;
        // auto p3cm = p3.Boost(-beta);
        // double p4Mag = p4.P();
        // double pot4 = Potential(p4Mag, rho);
        // double m4U2 = m4*m4+pow(pot4, 2)+2*pot4*sqrt(p4Mag*p4Mag+m4*m4);
        // return p1.E() + p2.E() - p3.E() - sqrt((p1 + p2 - p3).P2() + m4*m4);
        return s - 2*(p1+p2)*p3 + m3*m3 - m4*m4;
    };

    FourVector dummy(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta, 1);
    std::vector<double> guess;
    double a = (p1.E() + p2.E())/(2*m3);
    double b = -(p1.P()*p1.CosAngle(dummy)+p2.P()*p2.CosAngle(dummy));
    double c = m3*(p1.E() + p2.E())-s/2;
    double det = b*b-4*a*c;
    double rel_guess = s/(2*(p1.E()+p2.E()-(p1.P()*p1.CosAngle(dummy)+p2.P()*p2.CosAngle(dummy))));
    guess.emplace_back(det > 0 ? (-b + sqrt(det))/(2*a) : 0);
    guess.emplace_back(rel_guess);
    spdlog::info("Guess = {}", fmt::join(guess, ", ")); 
    spdlog::info("Func(guess) = {}, {}", _func(guess[0]), _func(guess[1]));
    double x = guess[1];
    double result_old = -std::numeric_limits<double>::infinity();
    while(x > -guess[1] && _func(guess[0])*_func(guess[1]) > 0) {
        double result = _func(x);
        spdlog::info("Func({}) = {}", x, _func(x));
        x -= 10;
        if(result * _func(guess[1]) < 0) break;
        if(result < result_old) break;
        result_old = result;
    }

    Brent brent(_func);
    double p3Mag = brent.CalcRoot(guess[0], guess[1]);
    spdlog::info("p3mag = {}", p3Mag);
    double pot3 = 0; // Potential(p3Mag, rho);
    double energy = sqrt(p3Mag*p3Mag + m3*m3) + pot3;
    return {p3Mag*sinTheta*cos(phi), p3Mag*sinTheta*sin(phi), p3Mag*cosTheta, energy};
}
