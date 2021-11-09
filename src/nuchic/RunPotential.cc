#include "nuchic/RunPotential.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Potential.hh"
#include "nuchic/SymplecticIntegrator.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Random.hh"
#include "nuchic/Utilities.hh"

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <chrono>
#include <fstream>

double Hamiltonian(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p,
                   std::shared_ptr<nuchic::Potential> potential) {
    auto vals = potential -> operator()(p.P(), q.P());
    auto mass_eff = nuchic::Constant::mN + vals.rscalar + std::complex<double>(0, 1)*vals.iscalar;
    return sqrt(p.P2() + pow(mass_eff, 2)).real() + vals.rvector;
}

nuchic::ThreeVector dHamiltonian_dp(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p,
                                    std::shared_ptr<nuchic::Potential> potential) {
    auto vals = potential -> operator()(p.P(), q.P());
    auto dpot_dp = potential -> derivative_p(p.P(), q.P(), 0.1*p.P());

    auto mass_eff = nuchic::Constant::mN + vals.rscalar + std::complex<double>(0, 1)*vals.iscalar;
    double numerator = (vals.rscalar + nuchic::Constant::mN)*dpot_dp.rscalar + p.P();
    double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
    return numerator/denominator * p/p.P() + dpot_dp.rvector * p/p.P();
}

nuchic::ThreeVector dHamiltonian_dr(const nuchic::ThreeVector &q, const nuchic::ThreeVector &p,
                                    std::shared_ptr<nuchic::Potential> potential) {
    auto vals = potential -> operator()(p.P(), q.P());
    auto dpot_dr = potential -> derivative_r(p.P(), q.P(), 0.1*q.P());

    auto mass_eff = nuchic::Constant::mN + vals.rscalar + std::complex<double>(0, 1)*vals.iscalar;
    double numerator = (vals.rscalar + nuchic::Constant::mN)*dpot_dr.rscalar;
    double denominator = sqrt(pow(mass_eff, 2) + p.P2()).real();
    return numerator/denominator * q/q.P() + dpot_dr.rvector * q/q.P();
}

void RunPropagation(std::shared_ptr<nuchic::Potential> potential,
                                    std::shared_ptr<nuchic::Nucleus> nucleus,
                                    const YAML::Node &config) {
    auto kick_mom = config["KickMomentum"].as<std::vector<double>>();
    auto nevents = config["NEvents"].as<size_t>();

    double current_mom = kick_mom[0];
    constexpr double omega = 20;
    constexpr size_t time_steps = 10000;
    constexpr double step_size = 0.01;
    double r0 = config["r0"].as<double>();

    auto mom = potential -> BindingMomentum(r0);

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
            double costheta = nuchic::Random::Instance().Uniform(-1.0, 1.0);
            double sintheta = sqrt(1-costheta*costheta);
            double phi = nuchic::Random::Instance().Uniform(0.0, 2*M_PI);
            nuchic::ThreeVector q{r0, 0, 0};
            nuchic::ThreeVector p{current_mom*sintheta*cos(phi), current_mom*sintheta*sin(phi), current_mom*costheta};
            nuchic::SymplecticIntegrator si(q, p, potential, dHamiltonian_dr, dHamiltonian_dp, omega); 
            for(size_t j = 0; j < time_steps; ++j) {
                si.Step<2>(step_size);
                if(si.Q().Magnitude() > 6.0) {
                    escaped++;
                    average_steps += j;
                    break;
                }
            }
        }
        auto result = static_cast<double>(escaped)/static_cast<double>(nevents);
        auto avg_steps = static_cast<double>(average_steps)/static_cast<double>(escaped);
        fmt::print("Escaped/Total = {}, Average Steps to escape = {}\n", result, avg_steps);
        out << fmt::format("{:8.3f},{:8.3f},{:8.3f}\n", current_mom, result, avg_steps); 
        current_mom += kick_mom[2];
    }
    out.close();
}

void RunBinding(std::shared_ptr<nuchic::Potential> potential, const YAML::Node &config) {
    fmt::print("Potential binding\n");
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream out(filename);
    auto radii = config["Radii"].as<std::vector<double>>();
    double current_radius = radii[0];
    while(current_radius <= radii[1]) {
        double mom = potential -> BindingMomentum(current_radius);
        fmt::print("  Radius = {}, Binding = {}\n", current_radius, mom);
        out << fmt::format("{:8.3f},{:8.3f}\n", current_radius, mom);
        current_radius += radii[2];
    }
    out.close();
}


void RunEnergySpectrum(std::shared_ptr<nuchic::Potential> potential, const YAML::Node &config) {
    fmt::print("Potential energy spectrum\n");
    std::string filename = fmt::format("{}.dat", config["SaveAs"].as<std::string>());
    std::ofstream out(filename);
    auto radii = config["Radii"].as<std::vector<double>>();
    double current_radius = radii[0];
    auto kick_mom = config["KickMomentum"].as<std::vector<double>>();

    while(current_radius <= radii[1]) {
	 double current_mom = kick_mom[0];
	 while(current_mom <= kick_mom[1]) {
            double energy_free= pow(current_mom,2)/2/nuchic::Constant::mN;
            if(config["Potential"].as<std::string>() == "Cooper")
		     energy_free= sqrt(pow(current_mom,2)+pow(nuchic::Constant::mN,2))-nuchic::Constant::mN;
            double energy = potential -> EnergySpectrum (current_radius, current_mom);
            fmt::print("  Radius = {}, Momentum = {}, Energy = {}, Free Energy = {}\n", current_radius, current_mom, energy, energy_free);
            out << fmt::format("{:8.3f},{:8.3f},{:8.3f},{:8.3f}\n", current_radius, current_mom, energy, energy_free);
            current_mom += kick_mom[2];
	 }
        current_radius += radii[2];
    }
    out.close();
}



void nuchic::RunPotential(const std::string &runcard) {
    auto config = YAML::LoadFile(runcard);
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["seed"])
        seed = config["Initialize"]["seed"].as<unsigned int>();
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
        potential = std::make_shared<SchroedingerPotential>(nucleus, config["Schrodinger"].as<size_t>());

    // Generate events
    if(config["Mode"].as<std::string>() == "Propagation")
        RunPropagation(potential, nucleus, config);
    else if(config["Mode"].as<std::string>() == "Binding")
        RunBinding(potential, config);
    else if(config["Mode"].as<std::string>() == "EnergySpectrum")
        RunEnergySpectrum(potential, config);
}
