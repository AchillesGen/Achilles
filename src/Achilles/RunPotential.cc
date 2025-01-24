#include "Achilles/RunPotential.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Random.hh"
#include "Achilles/SymplecticIntegrator.hh"
#include "Achilles/Utilities.hh"

#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <chrono>
#include <fstream>

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
    auto dpot_dp = potential->derivative_p(p.P(), q.P(), 0.1 * p.P());

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
    auto dpot_dr = potential->derivative_r(p.P(), q.P(), 0.1 * q.P());

    auto mass_eff =
        achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
    double numerator = (vals.rscalar + achilles::Constant::mN) * dpot_dr.rscalar;
    double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
    return numerator / denominator * q / q.P() + dpot_dr.rvector * q / q.P();
}

void RunPropagation(std::shared_ptr<achilles::Potential> potential,
                    std::shared_ptr<achilles::Nucleus> nucleus, const YAML::Node &config) {
    auto kick_mom = config["KickMomentum"].as<std::vector<double>>();
    auto nevents = config["NEvents"].as<size_t>();

    double current_mom = kick_mom[0];
    constexpr double omega = 20;
    constexpr size_t time_steps = 10000;
    constexpr double step_size = 0.01;
    double r0 = config["r0"].as<double>();

    auto mom = potential->BindingMomentum(r0);

    fmt::print("Potential propagation running with r0={}\n", r0);
    fmt::print("  Generating {} events per momentum point\n", nevents);
    fmt::print("  Fermi momentum = {}\n", nucleus->FermiMomentum(r0));
    fmt::print("  Minimum binding momentum = {}\n", mom);
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream out(filename);
    while(current_mom <= kick_mom[1]) {
        fmt::print("  Kick momentum: {} MeV    ", current_mom);
        size_t escaped = 0;
        size_t average_steps = 0;
        for(size_t i = 0; i < nevents; ++i) {
            double costheta = achilles::Random::Instance().Uniform(-1.0, 1.0);
            double sintheta = sqrt(1 - costheta * costheta);
            double phi = achilles::Random::Instance().Uniform(0.0, 2 * M_PI);
            achilles::ThreeVector q{r0, 0, 0};
            achilles::ThreeVector p{current_mom * sintheta * cos(phi),
                                    current_mom * sintheta * sin(phi), current_mom * costheta};
            achilles::SymplecticIntegrator si(q, p, potential, dHamiltonian_dr, dHamiltonian_dp,
                                              omega);
            for(size_t j = 0; j < time_steps; ++j) {
                si.Step<2>(step_size);
                if(si.Q().Magnitude() > 6.0) {
                    escaped++;
                    average_steps += j;
                    break;
                }
            }
        }
        auto result = static_cast<double>(escaped) / static_cast<double>(nevents);
        auto avg_steps = static_cast<double>(average_steps) / static_cast<double>(escaped);
        fmt::print("Escaped/Total = {}, Average Steps to escape = {}\n", result, avg_steps);
        out << fmt::format("{:8.3f},{:8.3f},{:8.3f}\n", current_mom, result, avg_steps);
        current_mom += kick_mom[2];
    }
    out.close();
}

void RunBinding(std::shared_ptr<achilles::Potential> potential, const YAML::Node &config) {
    fmt::print("Potential binding\n");
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream out(filename);
    auto radii = config["Radii"].as<std::vector<double>>();
    double current_radius = radii[0];
    while(current_radius <= radii[1]) {
        double mom = potential->BindingMomentum(current_radius);
        fmt::print("  Radius = {}, Binding = {}\n", current_radius, mom);
        out << fmt::format("{:8.3f},{:8.3f}\n", current_radius, mom);
        current_radius += radii[2];
    }
    out.close();
}

void RunEnergySpectrum(std::shared_ptr<achilles::Potential> potential, const YAML::Node &config) {
    fmt::print("Potential energy spectrum\n");
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream out(filename);
    auto radii = config["Radii"].as<std::vector<double>>();
    double current_radius = radii[0];
    auto kick_mom = config["KickMomentum"].as<std::vector<double>>();

    while(current_radius <= radii[1]) {
        double current_mom = kick_mom[0];
        while(current_mom <= kick_mom[1]) {
            double energy_free = pow(current_mom, 2) / 2 / achilles::Constant::mN;
            if(config["Potential"].as<std::string>() == "Cooper")
                energy_free = sqrt(pow(current_mom, 2) + pow(achilles::Constant::mN, 2)) -
                              achilles::Constant::mN;
            double energy = potential->EnergySpectrum(current_radius, current_mom);
            fmt::print("  Radius = {}, Momentum = {}, Energy = {}, Free Energy = {}\n",
                       current_radius, current_mom, energy, energy_free);
            out << fmt::format("{:8.3f},{:8.3f},{:8.3f},{:8.3f}\n", current_radius, current_mom,
                               energy, energy_free);
            current_mom += kick_mom[2];
        }
        current_radius += radii[2];
    }
    out.close();
}

double binomial(size_t n, size_t k) {
    if(k > n)
        throw std::runtime_error(fmt::format("Binomial coefficient requires n > k: {}, {}", n, k));

    double result = 1.0;
    if(k > n - k) k = n - k;
    for(size_t i = 0; i < k; ++i) {
        result *= static_cast<double>(n - i);
        result /= static_cast<double>(i + 1);
    }

    return result;
}

double OscSeries(const double &eps, const std::function<double(size_t)> &func) {
    int steps = static_cast<int>(std::ceil(-1.31 * log10(eps)));
    fmt::print("NSteps = {}\n", steps);
    double d = pow(3.0 + sqrt(8), steps);
    d = (d + 1 / d) / 2;
    double b = -1, c = -d, s = 0;
    for(int i = 0; i < steps; ++i) {
        size_t idx = static_cast<size_t>(i);
        c = b - c;
        s += c * func(idx);
        b = (i + steps) * (i - steps) * b / ((i + 0.5) * (i + 1));
    }

    return s / d;
}

void RunApprox(std::shared_ptr<achilles::Potential> potential, const YAML::Node &config) {
    auto kick_mom = config["KickMomentum"].as<std::vector<double>>();
    auto current_mom = kick_mom[0];
    double r = config["r0"].as<double>();

    fmt::print("Potential Approximation running with r={}\n", r);
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream out(filename);
    while(current_mom <= kick_mom[1]) {
        auto p = current_mom;
        fmt::print("  Kick momentum: {} MeV    ", current_mom);
        double exact = potential->Hamiltonian(p, r);
        auto vals = potential->operator()(p, r);
        fmt::print("Exact = {}\n", exact);
        auto mass_eff =
            achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
        double approx = mass_eff.real() + vals.rvector;
        double approx2 = p + vals.rvector;
        fmt::print("Approx[-1] = {}\n", approx);
        auto func = [&](size_t k) {
            return (2 / (static_cast<double>(k) + 1) * binomial(2 * k, k) *
                    pow(p * p / (4.0 * mass_eff * mass_eff), static_cast<int>(k + 1)) * mass_eff)
                .real();
        };
        auto func2 = [&](size_t k) {
            return (2 / (static_cast<double>(k) + 1) * binomial(2 * k, k) *
                    pow(mass_eff * mass_eff / (4.0 * p * p), static_cast<int>(k + 1)) * p)
                .real();
        };
        approx += OscSeries(5e-4, func);
        approx2 += OscSeries(5e-4, func2);
        fmt::print("ApproxSeries = {}, {}\n", approx, approx2);
        out << fmt::format("{:8.3f},{:8.3f},{:8.3f},{:8.3f}\n", current_mom, exact, approx,
                           approx2);

        current_mom += kick_mom[2];
    }
}

void achilles::RunPotential(const std::string &runcard) {
    auto config = YAML::LoadFile(runcard);
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["seed"]) seed = config["Initialize"]["seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load setup
    auto nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Initialize Potential
    std::shared_ptr<Potential> potential;
    if(config["Potential"].as<std::string>() == "Cooper")
        potential = std::make_shared<CooperPotential>(nucleus);
    else if(config["Potential"].as<std::string>() == "Wiringa")
        potential = std::make_shared<WiringaPotential>(nucleus);
    else if(config["Potential"].as<std::string>() == "Schrodinger")
        potential =
            std::make_shared<SchroedingerPotential>(nucleus, config["Schrodinger"].as<size_t>());

    // Generate events
    if(config["Mode"].as<std::string>() == "Propagation")
        RunPropagation(potential, nucleus, config);
    else if(config["Mode"].as<std::string>() == "Binding")
        RunBinding(potential, config);
    else if(config["Mode"].as<std::string>() == "EnergySpectrum")
        RunEnergySpectrum(potential, config);
    else if(config["Mode"].as<std::string>() == "Approx")
        RunApprox(potential, config);
}
