#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include <regex>
#pragma GCC diagnostic pop

#include "spdlog/spdlog.h"

#include "Achilles/Constants.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

const std::map<std::size_t, std::string> Nucleus::ZToName = {
    {0, "mfp"}, {1, "H"},   {2, "He"},  {3, "Li"},  {6, "C"},
    {8, "O"},   {13, "Al"}, {18, "Ar"}, {20, "Ca"}, {26, "Fe"},
};

Nucleus::Nucleus(const std::size_t &Z, const std::size_t &A, const double &bEnergy,
                 const double &kf, const std::string &densityFilename, const FermiGas &fg,
                 std::unique_ptr<Density> _density)
    : binding(bEnergy), fermiMomentum(kf), fermi_gas(fg), density(std::move(_density)) {
    Initialize(Z, A);
    // TODO: Refactor elsewhere in the code, maybe make dynamic?
    // spdlog::info("Nucleus: inferring nuclear radius using 0.16
    // nucleons/fm^3."); constexpr double nucDensity = 0.16; radius =
    // std::cbrt(static_cast<double>(A) / (4 / 3 * M_PI * nucDensity));

    std::ifstream densityFile(densityFilename);
    if(!densityFile.is_open())
        throw std::runtime_error(
            fmt::format("Nucleus: Density file {} does not exist.", densityFilename));
    std::string lineContent;

    constexpr size_t HeaderLength = 16;
    for(size_t i = 0; i < HeaderLength; ++i) { std::getline(densityFile, lineContent); }

    double radius_{}, density_{}, densityErr{};
    std::vector<double> vecRadius, vecDensity;
    constexpr double minDensity = 1E-6;
    while(densityFile >> radius_ >> density_ >> densityErr) {
        if(density_ < minDensity && radius == 0) radius = radius_;
        vecRadius.push_back(std::move(radius_));
        vecDensity.push_back(std::move(density_));
    }

    rhoInterp.SetData(vecRadius, vecDensity);
    rhoInterp.CubicSpline();

    // Ensure the number of protons and neutrons are correct
    // NOTE: This only is checked at startup, so if density returns a varying
    // number of nucleons it will not necessarily be caught
    auto particles = density->GetConfiguration();
    if(particles.size() != nnucleons)
        throw std::runtime_error("Invalid density function! Incorrect number of nucleons.");

    std::size_t nProtons = 0, nNeutrons = 0;
    for(auto particle : particles) {
        if(particle.ID() == PID::proton()) nProtons++;
        if(particle.ID() == PID::neutron()) nNeutrons++;
    }

    if(nProtons != NProtons() || nNeutrons != NNeutrons())
        throw std::runtime_error("Invalid density function! Incorrect number "
                                 "of protons or neutrons.");
}

void Nucleus::Initialize(size_t Z, size_t A) {
    if(Z > A) {
        std::string errorMsg = "Requires the number of protons to be less than the total";
        errorMsg += " number of nucleons. Got " + std::to_string(Z);
        errorMsg += " protons and " + std::to_string(A) + " nucleons";
        throw std::runtime_error(errorMsg);
    }

    nnucleons = A;
    nprotons = Z;
    nneutrons = A - Z;

    // Based on PDG Monte-Carlo PIDs
    // Nuclear codes are given as a 10 digit number:
    // +/- 10LZZZAAAI
    // L: number of strange baryons
    // Z: number of protons
    // A: number of nucleons
    // I: excited state (0 is ground state)
    static constexpr int IDBase = 1000000000;
    static constexpr int ZBase = 10000;
    static constexpr int ABase = 10;
    int ID = IDBase + ZBase * static_cast<int>(Z) + ABase * static_cast<int>(A);
    m_pid = PID{ID};

    // Set flag to handle special case of hydrogen
    if(Z == 1 && A == 1) is_hydrogen = true;
}

Particles Nucleus::GenerateConfig() {
    // Handle special case of hydrogen
    if(is_hydrogen) { return {{PID::proton(), {ParticleInfo(PID::proton()).Mass(), 0, 0, 0}}}; }

    // Get a configuration from the density function
    Particles particles = density->GetConfiguration();

    for(Particle &particle : particles) {
        // Set momentum for each nucleon
        auto mom3 = GenerateMomentum(particle.Position().Magnitude());
        double energy2 = pow(particle.Info().Mass(), 2); // Constant::mN*Constant::mN;
        for(auto mom : mom3) energy2 += mom * mom;
        particle.Momentum() = FourVector(sqrt(energy2), mom3[0], mom3[1], mom3[2]);

        // Ensure status is set to background
        particle.Status() = ParticleStatus::background;
    }
    return particles;
}

const std::array<double, 3> Nucleus::GenerateMomentum(const double &position) noexcept {
    std::array<double, 3> momentum{};
    //    momentum[0] = Random::Instance().Uniform(0.0, FermiMomentum(position)); //To sample on a
    //    sphere, need to take a cube-root.
    momentum[0] = FermiMomentum(position) * pow(Random::Instance().Uniform(0.0, 1.0), 1. / 3.);
    momentum[1] = std::acos(Random::Instance().Uniform(-1.0, 1.0));
    momentum[2] = Random::Instance().Uniform(0.0, 2 * M_PI);

    return ToCartesian(momentum);
}

Nucleus Nucleus::MakeNucleus(const std::string &name, const double &bEnergy,
                             const double &fermiMomentum, const std::string &densityFilename,
                             const FermiGas &fg, std::unique_ptr<Density> density) {
    const std::regex regex("([0-9]+)([a-zA-Z]+)");
    std::smatch match;

    if(std::regex_match(name, match, regex)) {
        const std::size_t nucleons = std::stoul(match[1].str());
        const std::size_t protons = NameToZ(match[2].str());
        spdlog::info("Nucleus: parsing nuclear name '{0}', expecting a density "
                     "with A={1} total nucleons and Z={2} protons.",
                     name, nucleons, protons);
        return Nucleus(protons, nucleons, bEnergy, fermiMomentum, densityFilename, fg,
                       std::move(density));
    }

    throw std::runtime_error(fmt::format("Invalid nucleus {}.", name));
}

std::size_t Nucleus::NameToZ(const std::string &name) {
    auto it =
        std::find_if(ZToName.begin(), ZToName.end(),
                     [&name](const std::pair<int, std::string> &p) { return p.second == name; });
    if(it == ZToName.end())
        throw std::runtime_error(fmt::format("Invalid nucleus: {} does not exist.", name));
    return it->first;
}

const std::string Nucleus::ToString() const noexcept {
    return std::to_string(NNucleons()) + ZToName.at(NProtons());
}

double Nucleus::FermiMomentum(const double &position) const noexcept {
    double rho = Rho(position);
    double result{};
    switch(fermi_gas.type) {
    // TODO: This doesn't give 20% of the nucleons with momentum greater than the Fermi momentum
    // It just changes the upper limit 20% of the time. Is this acceptable? Or do we need to change
    // how the momentum is generated?
    case FermiGasType::CorrelatedLocal: {
        double kF = std::cbrt(rho * 3 * M_PI * M_PI) * Constant::HBARC;
        if(Random::Instance().Uniform(0.0, 1.0) < fermi_gas.params[0]) {
            result = fermi_gas.params[1];
        } else {
            result = kF;
        }
    }
    case FermiGasType::Local:
        result = std::cbrt(rho * 3 * M_PI * M_PI) * Constant::HBARC;
        break;
    case FermiGasType::Global:
        //        static constexpr double small = 1E-2; //This leads to a large small momentum
        //        contribution at large radius, which is unrealistic, make small very snall ? result
        //        = rho < small ? small : fermiMomentum;
        result = fermiMomentum; // Better, but sampling is still srong, the pdf has p2 factor, need
                                // to take cube-root somewhere
        break;
    }

    return result;
}
