// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "catch2/catch_approx.hpp"
#include "catch2/catch_template_test_macros.hpp"
#include "catch2/catch_test_macros.hpp"

#include "Achilles/Constants.hh"
#include "Achilles/Potential.hh"
#include "Achilles/SymplecticIntegrator.hh"
#include "Achilles/ThreeVector.hh"
#include "mock_classes.hh"

#include <complex>

using achilles::ThreeVector;

double Rho(double r) {
    static constexpr double a = 0.09320982;
    static constexpr double b = 0.10665348;
    static constexpr double c = 0.42438208;
    static constexpr double d = 2.80880048;

    return a * pow(r, c) * exp(-b * pow(r, d));
}

double dRho(double r) {
    static constexpr double a = 0.09320982;
    static constexpr double b = 0.10665348;
    static constexpr double c = 0.42438208;
    static constexpr double d = 2.80880048;

    return a * pow(r, c - 1) * exp(-b * pow(r, d)) * (c - b * d * pow(r, d));
}

double Potential(double p, double r, double rho0) {
    const double rho_ratio = Rho(r) / rho0;
    const double alpha = 15.52 * rho_ratio + 24.93 * pow(rho_ratio, 2);
    const double beta = -116 * rho_ratio;
    const double lambda = (3.29 - 0.373 * rho_ratio) * achilles::Constant::HBARC;

    return alpha + beta / (1 + pow(p / lambda, 2));
}

double dPotential_dp(double p, double r, double rho0) {
    const double rho_ratio = Rho(r) / rho0;
    const double beta = -116 * rho_ratio;
    const double lambda = (3.29 - 0.373 * rho_ratio) * achilles::Constant::HBARC;

    return -2 * beta * p * pow(lambda, 2) / pow(p * p + pow(lambda, 2), 2);
}

double dPotential_dr(double p, double r, double rho0) {
    const double alpha0 = 15.52 / rho0;
    const double alpha1 = 24.93 / pow(rho0, 2);

    const double beta0 = -116 / rho0;
    const double lambda0 = 3.29 * achilles::Constant::HBARC;
    const double lambda1 = -0.373 / rho0 * achilles::Constant::HBARC;

    const double rho = Rho(r);

    const double term1 = alpha0;
    const double term2 = 2 * alpha1 * rho;
    const double term3 = beta0 / (1 + pow(p / (lambda0 - lambda1 * rho), 2));
    const double term4 =
        2 * beta0 * p * p * lambda1 * rho /
        (pow(lambda0 - lambda1 * rho, 3) * (1 + pow(p / (lambda0 - lambda1 * rho), 2)));

    return (term1 + term2 + term3 - term4) * dRho(r);
}

double Hamiltonian(const achilles::ThreeVector &q, const achilles::ThreeVector &p,
                   std::shared_ptr<achilles::Potential> potential) {
    auto vals = potential->operator()(p.P(), q.P());
    auto mass_eff =
        achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
    return sqrt(p.P2() + pow(mass_eff, 2)).real() + vals.rvector;
}

achilles::ThreeVector dHamiltonian_dp(const achilles::ThreeVector &q,
                                      const achilles::ThreeVector &p,
                                      std::shared_ptr<achilles::Potential> potential) {
    auto vals = potential->operator()(p.P(), q.P());
    auto dpot_dp = potential->derivative_p(p.P(), q.P());

    auto mass_eff =
        achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
    double numerator = (vals.rscalar + achilles::Constant::mN) * dpot_dp.rscalar + p.P();
    double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
    return numerator / denominator * p / p.P() + dpot_dp.rvector * p / p.P();
}

achilles::ThreeVector dHamiltonian_dr(const achilles::ThreeVector &q,
                                      const achilles::ThreeVector &p,
                                      std::shared_ptr<achilles::Potential> potential) {
    auto vals = potential->operator()(p.P(), q.P());
    auto dpot_dr = potential->derivative_r(p.P(), q.P());

    auto mass_eff =
        achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
    double numerator = (vals.rscalar + achilles::Constant::mN) * dpot_dr.rscalar;
    double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
    return numerator / denominator * q / q.P() + dpot_dr.rvector * q / q.P();
}

template <typename T> std::shared_ptr<T> MakePotential(std::shared_ptr<achilles::Nucleus> nuc) {
    return std::make_shared<T>(nuc);
}

TEMPLATE_TEST_CASE("Symplectic Integrator", "[Symplectic]", achilles::CooperPotential,
                   achilles::WiringaPotential) {
    constexpr double r0 = -1;
    constexpr double pmag = 275;
    achilles::ThreeVector q{r0, 0, 0};
    achilles::ThreeVector p{0, pmag, 0};
    constexpr size_t nsteps = 10000;
    constexpr double step_size = 0.01;
    constexpr double omega = 20;
    constexpr size_t AA = 12;

    auto nucleus = MakeNucleus(6, AA, 12, Rho);
    auto potential = MakePotential<TestType>(nucleus);

    auto dHamiltonian_dr_func = [&](const ThreeVector &p_, const ThreeVector &q_,
                                    std::shared_ptr<achilles::Potential> pot_) -> ThreeVector {
        return dHamiltonian_dr(p_, q_, pot_);
    };

    auto dHamiltonian_dp_func = [&](const ThreeVector &p_, const ThreeVector &q_,
                                    std::shared_ptr<achilles::Potential> pot_) -> ThreeVector {
        return dHamiltonian_dp(p_, q_, pot_);
    };

    const double E0 = Hamiltonian(q, p, potential);

    // The defining property of a symplectic integrator is that it conserves the
    // Hamiltonian (the total energy) to within a bounded error that shrinks as the
    // integration order increases. We integrate the trajectory and check that the
    // relative energy drift stays small, rather than dumping the trajectory to a
    // file for visual inspection.
    SECTION("Order 2 conserves energy") {
        achilles::SymplecticIntegrator si(q, p, potential, dHamiltonian_dr_func,
                                          dHamiltonian_dp_func, omega);
        for(size_t i = 0; i < nsteps; ++i) si.Step<2>(step_size);

        const double Ef = Hamiltonian(si.Q(), si.P(), potential);
        CHECK(std::abs(Ef - E0) / E0 < 1e-3);

        // The auxiliary (x, y) copy tracked by the integrator must stay locked to
        // the primary (q, p) phase-space point.
        CHECK(Hamiltonian(si.State().x, si.State().y, potential) ==
              Catch::Approx(Ef).epsilon(1e-6));
    }

    SECTION("Order 4 conserves energy better than order 2") {
        achilles::SymplecticIntegrator si2(q, p, potential, dHamiltonian_dr_func,
                                           dHamiltonian_dp_func, omega);
        achilles::SymplecticIntegrator si4(q, p, potential, dHamiltonian_dr_func,
                                           dHamiltonian_dp_func, omega);
        for(size_t i = 0; i < nsteps; ++i) {
            si2.Step<2>(step_size);
            si4.Step<4>(step_size);
        }

        const double drift2 = std::abs(Hamiltonian(si2.Q(), si2.P(), potential) - E0) / E0;
        const double drift4 = std::abs(Hamiltonian(si4.Q(), si4.P(), potential) - E0) / E0;
        CHECK(drift4 < 1e-4);
        CHECK(drift4 < drift2);
    }

    SECTION("Order 6 conserves energy better than order 2") {
        // Orders 4 and 6 both reach the ~1e-7 error floor for these parameters, so we
        // do not compare them against each other; instead we check order 6 against the
        // much larger order-2 drift, which is a stable, well-separated comparison.
        achilles::SymplecticIntegrator si2(q, p, potential, dHamiltonian_dr_func,
                                           dHamiltonian_dp_func, omega);
        achilles::SymplecticIntegrator si6(q, p, potential, dHamiltonian_dr_func,
                                           dHamiltonian_dp_func, omega);
        for(size_t i = 0; i < nsteps; ++i) {
            si2.Step<2>(step_size);
            si6.Step<6>(step_size);
        }

        const double drift2 = std::abs(Hamiltonian(si2.Q(), si2.P(), potential) - E0) / E0;
        const double drift6 = std::abs(Hamiltonian(si6.Q(), si6.P(), potential) - E0) / E0;
        CHECK(drift6 < 1e-5);
        CHECK(drift6 < drift2);
    }
}
