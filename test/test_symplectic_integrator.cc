#include "catch2/catch.hpp"

#include "mock_classes.hh"
#include "nuchic/Constants.hh"
#include "nuchic/SymplecticIntegrator.hh"
#include "spdlog/spdlog.h"
#include "nuchic/Potential.hh"

#include <fstream>
#include <complex>

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

double Hamiltonian(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p,
                   std::shared_ptr<nuchic::Potential> potential) {
    auto vals = potential -> operator()(p.P(), q.P());
    auto mass_eff = nuchic::Constant::mN + vals.rscalar + std::complex<double>(0, 1)*vals.iscalar;
    return sqrt(p.P2() + pow(mass_eff, 2)).real() + vals.rvector;
}

nuchic::ThreeVector dHamiltonian_dp(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p,
                                    std::shared_ptr<nuchic::Potential> potential) {
    auto vals = potential -> operator()(p.P(), q.P());
    auto dpot_dp = potential -> derivative_p(p.P(), q.P());

    auto mass_eff = nuchic::Constant::mN + vals.rscalar + std::complex<double>(0, 1)*vals.iscalar;
    double numerator = (vals.rscalar + nuchic::Constant::mN)*dpot_dp.rscalar + p.P();
    double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
    return numerator/denominator * p/p.P() + dpot_dp.rvector * p/p.P();
}

nuchic::ThreeVector dHamiltonian_dr(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p,
                                    std::shared_ptr<nuchic::Potential> potential) {
    auto vals = potential -> operator()(p.P(), q.P());
    auto dpot_dr = potential -> derivative_r(p.P(), q.P());

    auto mass_eff = nuchic::Constant::mN + vals.rscalar + std::complex<double>(0, 1)*vals.iscalar;
    double numerator = (vals.rscalar + nuchic::Constant::mN)*dpot_dr.rscalar;
    double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
    return numerator/denominator * q/q.P() + dpot_dr.rvector * q/q.P();
}

template<typename T>
std::shared_ptr<T> MakePotential(std::shared_ptr<nuchic::Nucleus> nuc) {
    return std::make_shared<T>(nuc);
}

TEMPLATE_TEST_CASE("Symplectic Integrator", "[Symplectic]", nuchic::CooperPotential, nuchic::WiringaPotential) {
    constexpr double r0 = -1;
    constexpr double pmag = 275;
    nuchic::ThreeVector q{r0, 0, 0};
    nuchic::ThreeVector p{0, pmag, 0};
    constexpr size_t nsteps = 10000;
    constexpr double step_size = 0.01;
    constexpr double omega = 20;
    constexpr size_t AA = 12;

    auto nucleus = std::make_shared<MockNucleus>();
    ALLOW_CALL(*nucleus, NNucleons())
        .LR_RETURN((AA));
    ALLOW_CALL(*nucleus, Rho(trompeloeil::gt(0)))
        .LR_RETURN((Rho(_1)));
    auto potential = MakePotential<TestType>(nucleus);

    nuchic::SymplecticIntegrator si(q, p, potential, dHamiltonian_dr, dHamiltonian_dp, omega);

    SECTION("Order 2") {
        std::ofstream out("symplectic2_" + potential -> Name() + ".txt");
        const double E0 = Hamiltonian(q, p, potential);
        spdlog::info("Initial Hamiltonian Value: {}", Hamiltonian(q, p, potential));
        out << "X,Y,Z,Px,Py,Pz,E\n";
        out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
        out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << E0 << "\n";
        for(size_t i = 0; i < nsteps; ++i) {
            si.Step<2>(step_size);
            const double Ei = Hamiltonian(si.Q(), si.P(), potential);
            out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
            out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << Ei << "\n";
        }
        const double Ef = Hamiltonian(si.Q(), si.P(), potential);
        const double Ef2 = Hamiltonian(si.State().x, si.State().y, potential);
        spdlog::info("Final Hamiltonian Value: {}, {}, {}, {}",
                Hamiltonian(si.Q(), si.P(), potential), std::abs(Ef-E0)/E0, Ef2, std::abs(Ef2-E0)/E0);
        spdlog::info("Final Position: {}, {}, {}", si.Q(), si.State().x, (si.Q() - si.State().x).P()/si.Q().P());
        spdlog::info("Final Momentum: {}, {}, {}", si.P(), si.State().y, (si.P() - si.State().y).P()/si.P().P());
        out.close();
    }

    SECTION("Order 4") {
        std::ofstream out("symplectic4_" + potential -> Name() + ".txt");
        const double E0 = Hamiltonian(q, p, potential);
        spdlog::info("Initial Hamiltonian Value: {}", Hamiltonian(q, p, potential));
        out << "X,Y,Z,Px,Py,Pz,E\n";
        out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
        out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << E0 << "\n";
        for(size_t i = 0; i < nsteps; ++i) {
            si.Step<4>(step_size);
            const double Ei = Hamiltonian(si.Q(), si.P(), potential);
            out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
            out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << Ei << "\n";
        }
        const double Ef = Hamiltonian(si.Q(), si.P(), potential);
        spdlog::info("Final Hamiltonian Value: {}, {}", Hamiltonian(si.Q(), si.P(), potential), std::abs(Ef-E0)/E0);
        spdlog::info("Final Position: {}, {}", si.Q(), si.State().x);
        spdlog::info("Final Momentum: {}, {}", si.P(), si.State().y);
        out.close();
    }

    SECTION("Order 6") {
        std::ofstream out("symplectic6_" + potential -> Name() + ".txt");
        const double E0 = Hamiltonian(q, p, potential);
        spdlog::info("Initial Hamiltonian Value: {}", Hamiltonian(q, p, potential));
        out << "X,Y,Z,Px,Py,Pz,E\n";
        out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
        out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << E0 << "\n";
        for(size_t i = 0; i < nsteps; ++i) {
            si.Step<6>(step_size);
            const double Ei = Hamiltonian(si.Q(), si.P(), potential);
            out << si.Q().X() << "," << si.Q().Y() << "," << si.Q().Z() << ",";
            out << si.P().X() << "," << si.P().Y() << "," << si.P().Z() << "," << Ei << "\n";
        }
        const double Ef = Hamiltonian(si.Q(), si.P(), potential);
        spdlog::info("Final Hamiltonian Value: {}, {}", Hamiltonian(si.Q(), si.P(), potential), std::abs(Ef-E0)/E0);
        spdlog::info("Final Position: {}, {}", si.Q(), si.State().x);
        spdlog::info("Final Momentum: {}, {}", si.P(), si.State().y);
        out.close();
    }
}
