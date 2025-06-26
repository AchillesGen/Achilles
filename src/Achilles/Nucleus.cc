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
    {0, "N"}, {1, "H"},   {2, "He"},  {3, "Li"},  {6, "C"},
    {8, "O"}, {13, "Al"}, {18, "Ar"}, {20, "Ca"}, {26, "Fe"},
};

Nucleus::Nucleus(const std::size_t &Z, const std::size_t &A, const double &bEnergy,
                 const double &kf, const std::string &protondensityFilename,
                 const std::string &neutrondensityFilename, const FermiGas &fg,
                 std::unique_ptr<Density> _density)
    : binding(bEnergy), fermiMomentum(kf), fermi_gas(fg), density(std::move(_density)) {
    Initialize(Z, A);
    // TODO: Refactor elsewhere in the code, maybe make dynamic?
    // spdlog::info("Nucleus: inferring nuclear radius using 0.16
    // nucleons/fm^3."); constexpr double nucDensity = 0.16; radius =
    // std::cbrt(static_cast<double>(A) / (4 / 3 * M_PI * nucDensity));

    // Density for protons
    auto protondensityPathFile = Filesystem::FindFile(protondensityFilename, "Nucleus");
    std::ifstream protondensityFile(protondensityPathFile);
    if(!protondensityFile.is_open())
        throw std::runtime_error(
            fmt::format("Nucleus: Issue opening file {}.", protondensityPathFile));
    std::string lineContent;

    constexpr size_t HeaderLength = 16;
    for(size_t i = 0; i < HeaderLength; ++i) { std::getline(protondensityFile, lineContent); }

    double protonradius_{}, protondensity_{}, protondensityErr{};
    std::vector<double> protonvecRadius, protonvecDensity;
    constexpr double minDensity = 1E-6;
    while(protondensityFile >> protonradius_ >> protondensity_ >> protondensityErr) {
        if(protondensity_ < minDensity && radius == 0) radius = protonradius_;
        protonvecRadius.push_back(std::move(protonradius_));
        protonvecDensity.push_back(std::move(protondensity_));
    }

    protonrhoInterp.SetData(protonvecRadius, protonvecDensity);
    protonrhoInterp.CubicSpline();

    // Density for neutrons
    auto neutrondensityPathFile = Filesystem::FindFile(neutrondensityFilename, "Nucleus");
    std::ifstream neutrondensityFile(neutrondensityPathFile);
    if(!neutrondensityFile.is_open())
        throw std::runtime_error(
            fmt::format("Nucleus: Issue opening file {}.", neutrondensityPathFile));
    // std::string lineContent;

    // constexpr size_t HeaderLength = 16;
    for(size_t i = 0; i < HeaderLength; ++i) { std::getline(neutrondensityFile, lineContent); }

    double neutronradius_{}, neutrondensity_{}, neutrondensityErr{};
    std::vector<double> neutronvecRadius, neutronvecDensity;
    // constexpr double minDensity = 1E-6;
    while(neutrondensityFile >> neutronradius_ >> neutrondensity_ >> neutrondensityErr) {
        // if(density_ < minDensity && radius == 0) radius = radius_;
        neutronvecRadius.push_back(std::move(neutronradius_));
        neutronvecDensity.push_back(std::move(neutrondensity_));
    }

    neutronrhoInterp.SetData(neutronvecRadius, neutronvecDensity);
    neutronrhoInterp.CubicSpline();

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

    // Set flag to handle special case of hydrogen or free neutron
    if(Z == 1 && A == 1) is_hydrogen = true;
    if(Z == 0 && A == 1) is_free_neutron = true;
}

Particles Nucleus::GenerateConfig() {
    // Handle special case of hydrogen
    if(is_hydrogen) {
        return {Particle{PID::proton(), {ParticleInfo(PID::proton()).Mass(), 0, 0, 0}}};
    }
    if(is_free_neutron) {
        return {Particle{PID::neutron(), {ParticleInfo(PID::neutron()).Mass(), 0, 0, 0}}};
    }

    // Get a configuration from the density function
    Particles particles = density->GetConfiguration();

    for(Particle &particle : particles) {
        // Set momentum for each nucleon
        auto mom3 = GenerateMomentum(particle.Position().Magnitude(),particle.ID());
        double energy2 = pow(particle.Info().Mass(), 2); // Constant::mN*Constant::mN;
        for(auto mom : mom3) energy2 += mom * mom;
        particle.Momentum() = FourVector(sqrt(energy2), mom3[0], mom3[1], mom3[2]);

        // Ensure status is set to background
        particle.Status() = ParticleStatus::background;
    }
    return particles;
}

const std::array<double, 3> Nucleus::GenerateMomentum(const double &position, const PID &pid) noexcept {
    std::array<double, 3> momentum{};
    momentum[0] = SampleMagnitudeMomentum(position,pid);
    momentum[1] = std::acos(Random::Instance().Uniform(-1.0, 1.0));
    momentum[2] = Random::Instance().Uniform(0.0, 2 * M_PI);

    return ToCartesian(momentum);
}

double Nucleus::SampleMagnitudeMomentum(const double &position, const PID &pid) noexcept {
    // NOTE: To sample on a sphere, need to take a cube-root.
    double kf = FermiMomentum(position,pid);
    if(fermi_gas.correlated) {
        if(Random::Instance().Uniform(0.0, 1.0) > fermi_gas.SRCfraction) {
            return kf * std::cbrt(Random::Instance().Uniform(0.0, 1.0));
        } else {
            double x = Random::Instance().Uniform(0.0, 1.0);
            return kf / (1. + 1. / fermi_gas.lambdaSRC - x);
        }
    }
    return kf * std::cbrt(Random::Instance().Uniform(0.0, 1.0));
}

Nucleus Nucleus::MakeNucleus(const std::string &name, const double &bEnergy,
                             const double &fermiMomentum, const std::string &protondensityFilename,
                             const std::string &neutrondensityFilename, const FermiGas &fg,
                             std::unique_ptr<Density> density) {
    const std::regex regex("([0-9]+)([a-zA-Z]+)");
    std::smatch match;

    if(std::regex_match(name, match, regex)) {
        const std::size_t nucleons = std::stoul(match[1].str());
        const std::size_t protons = NameToZ(match[2].str());
        spdlog::info("Nucleus: parsing nuclear name '{0}', expecting a density "
                     "with A={1} total nucleons and Z={2} protons.",
                     name, nucleons, protons);
        return Nucleus(protons, nucleons, bEnergy, fermiMomentum, protondensityFilename,
                       neutrondensityFilename, fg, std::move(density));
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

double Nucleus::FermiMomentum(const double &position, const PID &nuc_pid) const {
    double rho = 0.;
    if (nuc_pid == PID::proton()) rho = ProtonRho(position);
    else if (nuc_pid == PID::neutron()) rho = NeutronRho(position);
    else {
        throw std::runtime_error(fmt::format("Fermi Momentum for: {} does not exist.", nuc_pid));
    }

    //double rho = Rho(position);
    double result{};
    switch(fermi_gas.type) {
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
