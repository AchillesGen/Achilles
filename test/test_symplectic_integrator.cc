#include "catch2/catch.hpp"

#include "nuchic/Constants.hh"
#include "nuchic/SymplecticIntegrator.hh"
#include "spdlog/spdlog.h"

#include <fstream>

constexpr double a = 0.09320982;
constexpr double b = 0.10665348;
constexpr double c = 0.42438208;
constexpr double d = 2.80880048;

double Rho(double r) {
    return a*pow(r, c)*exp(-b*pow(r, d));
}

double dRho(double r) {
    return a*pow(r, c-1)*exp(-b*pow(r, d))*(c-b*d*pow(r, d));
}

double Potential(double p, double r, double rho0) {
    const double rho_ratio = Rho(r)/rho0;
    const double alpha = 15.52*rho_ratio + 24.93*pow(rho_ratio, 2);
    const double beta = -116*rho_ratio;
    const double lambda = (3.29 - 0.373*rho_ratio)*nuchic::Constant::HBARC;

    return alpha + beta/(1+pow(p/lambda, 2));
}

double dPotential_dp(double p, double r, double rho0) {
    const double rho_ratio = Rho(r)/rho0;
    const double beta = -116*rho_ratio;
    const double lambda = (3.29 - 0.373*rho_ratio)*nuchic::Constant::HBARC;

    return -2*beta*p*pow(lambda, 2)/pow(p*p+pow(lambda, 2), 2);
}

double dPotential_dr(double p, double r, double rho0) {
    const double alpha0 = 15.52/rho0;
    const double alpha1 = 24.93/pow(rho0, 2);

    const double beta0 = -116/rho0;
    const double lambda0 = 3.29*nuchic::Constant::HBARC;
    const double lambda1 = -0.373/rho0*nuchic::Constant::HBARC;

    const double rho = Rho(r);

    const double term1 = alpha0;
    const double term2 = 2*alpha1*rho;
    const double term3 = beta0/(1+pow(p/(lambda0-lambda1*rho), 2));
    const double term4 = 2*beta0*p*p*lambda1*rho/(pow(lambda0-lambda1*rho, 3)*(1+pow(p/(lambda0-lambda1*rho), 2)));

    return (term1+term2+term3-term4)*dRho(r);
}

double Hamiltonian(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p) {
    const double Tp = sqrt(p.P2() + pow(nuchic::Constant::mN, 2)) - nuchic::Constant::mN;
    return sqrt(p.P2() + pow(nuchic::Constant::mN, 2)) + Potential(Tp, q.P(), 0.16);
}

nuchic::ThreeVector dHamiltonian_dp(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p) {
    const double E = sqrt(p.P2() + pow(nuchic::Constant::mN, 2));
    const double Tp = sqrt(p.P2() + pow(nuchic::Constant::mN, 2)) - nuchic::Constant::mN;
    return p/E + dPotential_dp(Tp, q.P(), 0.16)*p/p.P();
}

nuchic::ThreeVector dHamiltonian_dr(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p) {
    const double Tp = sqrt(p.P2() + pow(nuchic::Constant::mN, 2)) - nuchic::Constant::mN;
    return dPotential_dr(Tp, q.P(), 0.16)*q/q.P();
}

TEST_CASE("Symplectic Integrator", "[Symplectic]") {
    constexpr double r0 = -1;
    constexpr double pmag = 300;
    nuchic::ThreeVector q{r0, 0, 0};
    nuchic::ThreeVector p{0, pmag, 0};
    constexpr size_t nsteps = 40000;
    constexpr double step_size = 0.01;
    constexpr double omega = 20;

    nuchic::SymplecticIntegrator si(q, p, dHamiltonian_dr, dHamiltonian_dp, omega);

    SECTION("Order 2") {
        std::ofstream out("symplectic2.txt");
        const double E0 = Hamiltonian(q, p);
        spdlog::info("Initial Hamiltonian Value: {}", Hamiltonian(q, p));
        out << "X,Y,Z,Px,Py,Pz,E\n";
        out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
        out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << E0 << "\n";
        for(size_t i = 0; i < nsteps; ++i) {
            si.Step<2>(step_size);
            const double Ei = Hamiltonian(si.Q(), si.P());
            out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
            out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << Ei << "\n";
        }
        const double Ef = Hamiltonian(si.Q(), si.P());
        const double Ef2 = Hamiltonian(si.State().x, si.State().y);
        spdlog::info("Final Hamiltonian Value: {}, {}, {}, {}",
                Hamiltonian(si.Q(), si.P()), std::abs(Ef-E0)/E0, Ef2, std::abs(Ef2-E0)/E0);
        spdlog::info("Final Position: {}, {}, {}", si.Q(), si.State().x, (si.Q() - si.State().x).P()/si.Q().P());
        spdlog::info("Final Momentum: {}, {}, {}", si.P(), si.State().y, (si.P() - si.State().y).P()/si.P().P());
        out.close();
    }

    SECTION("Order 4") {
        std::ofstream out("symplectic4.txt");
        const double E0 = Hamiltonian(q, p);
        spdlog::info("Initial Hamiltonian Value: {}", Hamiltonian(q, p));
        out << "X,Y,Z,Px,Py,Pz,E\n";
        out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
        out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << E0 << "\n";
        for(size_t i = 0; i < nsteps; ++i) {
            si.Step<4>(step_size);
            const double Ei = Hamiltonian(si.Q(), si.P());
            out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
            out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << Ei << "\n";
        }
        const double Ef = Hamiltonian(si.Q(), si.P());
        spdlog::info("Final Hamiltonian Value: {}, {}", Hamiltonian(si.Q(), si.P()), std::abs(Ef-E0)/E0);
        spdlog::info("Final Position: {}, {}", si.Q(), si.State().x);
        spdlog::info("Final Momentum: {}, {}", si.P(), si.State().y);
        out.close();
    }

    SECTION("Order 6") {
        std::ofstream out("symplectic6.txt");
        const double E0 = Hamiltonian(q, p);
        spdlog::info("Initial Hamiltonian Value: {}", Hamiltonian(q, p));
        out << "X,Y,Z,Px,Py,Pz,E\n";
        out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
        out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << E0 << "\n";
        for(size_t i = 0; i < nsteps; ++i) {
            si.Step<6>(step_size);
            const double Ei = Hamiltonian(si.Q(), si.P());
            out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
            out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << Ei << "\n";
        }
        const double Ef = Hamiltonian(si.Q(), si.P());
        spdlog::info("Final Hamiltonian Value: {}, {}", Hamiltonian(si.Q(), si.P()), std::abs(Ef-E0)/E0);
        spdlog::info("Final Position: {}, {}", si.Q(), si.State().x);
        spdlog::info("Final Momentum: {}, {}", si.P(), si.State().y);
        out.close();
    }
}
