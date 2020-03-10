#include <cmath>
#include <map>
#include <regex>

#include <iostream>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "nuchic/Constants.hh"
#include "nuchic/ThreeVector.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Utilities.hh"

using namespace nuchic;

const std::map<std::size_t, std::string> Nucleus::ZToName = {
    {0, "mfp"},
    {1, "H"},
    {2, "He"},
    {3, "Li"},
    {6, "C"},
    {8, "O"},
    {13, "Al"},
    {18, "Ar"},
    {20, "Ca"},
    {26, "Fe"},
};

Nucleus::Nucleus(const std::size_t& Z, const std::size_t& A, const double& bEnergy,
                         const double& kf, const std::function<Particles()>& _density) 
                        : binding(bEnergy), fermiMomentum(kf), density(_density) {

    if(Z > A) {
        std::string errorMsg = "Requires the number of protons to be less than the total";
        errorMsg += " number of nucleons. Got " + std::to_string(Z);
        errorMsg += " protons and " + std::to_string(A) + " nucleons";
        throw std::runtime_error(errorMsg);
    }
    
    nucleons.resize(A);
    protons.resize(Z);
    neutrons.resize(A-Z);
    spdlog::info("Nucleus: inferring nuclear radius using 0.16 nucleons/fm^3.");
    radius = pow(static_cast<double>(A) / (4.0 / 3.0 * M_PI * 0.16), 1.0 / 3.0);
    potential = sqrt(Constant::mN*Constant::mN 
                     + pow(fermiMomentum, 2)) - Constant::mN + 8;
}

void Nucleus::SetNucleons(Particles& _nucleons) noexcept {
    nucleons = _nucleons;
    std::size_t proton_idx = 0;
    std::size_t neutron_idx = 0;
    for(auto particle : nucleons) {
        if(particle.PID() == 2212) protons[proton_idx++] = particle;
        if(particle.PID() == 2112) neutrons[neutron_idx++] = particle;
    }
}

bool Nucleus::Escape(Particle& particle) noexcept {
    // Remove background particles
    if(particle.Status() == 0) return false;

    // Special case for testing pN cross-section
    if(particle.Status() == -2) return true;

    // Calculate kinetic energy, and if less than potential it is captured
    const double totalEnergy = sqrt(particle.Momentum().P2() + particle.Momentum().M2());
    const double kineticEnergy = totalEnergy - particle.Mass();
    if(kineticEnergy < potential) return false;

    // If the particle escapes, adjust momentum to account for this
    // TODO: This adjusts the mass. Is that acceptable?
    const double theta = particle.Momentum().Theta();
    const double phi = particle.Momentum().Phi();
    const double px = particle.Momentum().Px() - potential * std::sin(theta) * std::cos(phi);
    const double py = particle.Momentum().Py() - potential * std::sin(theta) * std::sin(phi);
    const double pz = particle.Momentum().Pz() - potential * std::cos(theta);
    particle.SetMomentum(FourVector(px, py, pz, particle.Momentum().E()));
    return true;
}

Particles Nucleus::GenerateConfig() {
    // Get a configuration from the density function
    Particles particles = density();

    // Ensure the number of protons and neutrons are correct
    if(particles.size() != nucleons.size())
        throw std::runtime_error("Invalid density function! Incorrect number of nucleons.");
    std::size_t nProtons = 0, nNeutrons = 0;
    for(Particle& particle : particles) {
        if(particle.PID() == 2212) nProtons++;
        if(particle.PID() == 2112) nNeutrons++;

        // Set momentum for each nucleon
        auto mom3 = GenerateMomentum();
        double energy2 = Constant::mN*Constant::mN;
        for(auto mom : mom3) energy2 += mom*mom;
        particle.SetMomentum(FourVector(mom3[0], mom3[1], mom3[2], sqrt(energy2)));
    }
    if(nProtons != NProtons() || nNeutrons != NNeutrons())
        throw std::runtime_error("Invalid density function! Incorrect number of protons and neutrons.");

    // Update the nucleons in the nucleus
    SetNucleons(particles);
    return particles;
}

const std::array<double, 3> Nucleus::GenerateMomentum() noexcept {
    std::array<double, 3> momentum;
    momentum[0] = rng.uniform(0.0, fermiMomentum);
    momentum[1] = std::acos(rng.uniform(-1.0, 1.0));
    momentum[2] = rng.uniform(0.0, 2*M_PI);

    return ToCartesian(momentum);
}

Nucleus Nucleus::MakeNucleus(const std::string& name, const double& bEnergy,
                                             const double& fermiMomentum,
                                             const std::function<Particles()>& density) {
    const std::regex regex("([0-9]+)([a-zA-Z]+)");
    std::smatch match;

    if(std::regex_match(name, match, regex)) {
        const std::size_t nucleons = std::stoul(match[1].str());
        const std::size_t protons = NameToZ(match[2].str()); 
        spdlog::info(
            "Nucleus: parsing nuclear name '{0}', expecting a density "
            "with A={1} total nucleons and Z={2} protons.", 
            name, nucleons, protons);
        return Nucleus(protons, nucleons, bEnergy, fermiMomentum, density);
    }

    throw std::runtime_error("Invalid nucleus " + name);
}

std::size_t Nucleus::NameToZ(const std::string& name) {
    auto it = std::find_if(ZToName.begin(), ZToName.end(),
                           [&name](const std::pair<int, std::string> &p) {
                               return p.second == name;
                          });
    if(it == ZToName.end()) 
        throw std::runtime_error("Invalid nucleus: " + name + " does not exist.");
    return it -> first;
}

const std::string Nucleus::ToString() const noexcept {
    return std::to_string(NNucleons()) + ZToName.at(NProtons());
}
