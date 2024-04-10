#include "Achilles/Event.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::Event;

Event::Event(Nucleus *nuc, std::vector<FourVector> mom, double vwgt)
    : m_nuc{nuc->NProtons(), nuc->NNeutrons()}, m_mom{std::move(mom)}, m_wgt{std::move(vwgt)} {
    m_hadrons = nuc->GenerateConfig();
}

Event::Event(const Event &other) {
    m_nuc = other.m_nuc;
    m_remnant = other.m_remnant;
    m_mom = other.m_mom;
    m_wgt = other.m_wgt;
    m_leptons = other.m_leptons;
    m_history = other.m_history;
    flux = other.flux;
}

void Event::Finalize() {
    size_t nA = 0, nZ = 0;
    for(const auto &part : m_hadrons) {
        if(part.IsExternal()) continue;
        if(part.ID() == PID::proton()) nZ++;
        nA++;
    }

    m_remnant = NuclearRemnant(nA, nZ);
}

void Event::Display() const {
    spdlog::trace("Leptons:");
    size_t idx = 0;
    for(const auto &particle : Leptons()) { spdlog::trace("\t{}: {}", ++idx, particle); }
    spdlog::trace("Hadrons:");
    idx = 0;
    for(const auto &particle : Hadrons()) { spdlog::trace("\t{}: {}", ++idx, particle); }
    spdlog::trace("Weight: {}", Weight());
}

achilles::vParticles Event::Particles() const {
    vParticles result;
    result.insert(result.end(), m_hadrons.begin(), m_hadrons.end());
    result.insert(result.end(), m_leptons.begin(), m_leptons.end());
    return result;
}

achilles::crefParticles Event::Protons(ParticleStatus status) const {
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::proton(); };
        return FilterParticles(m_hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::proton() && p.Status() == status;
    };
    return FilterParticles(m_hadrons, func);
}

achilles::crefParticles Event::Neutrons(ParticleStatus status) const {
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::neutron(); };
        return FilterParticles(m_hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::neutron() && p.Status() == status;
    };
    return FilterParticles(m_hadrons, func);
}

achilles::refParticles Event::Protons(ParticleStatus status) {
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::proton(); };
        return FilterParticles(m_hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::proton() && p.Status() == status;
    };
    return FilterParticles(m_hadrons, func);
}

achilles::refParticles Event::Neutrons(ParticleStatus status) {
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::neutron(); };
        return FilterParticles(m_hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::neutron() && p.Status() == status;
    };
    return FilterParticles(m_hadrons, func);
}

void Event::Rotate(const std::array<double, 9> &rot_mat) {
    for(auto &particle : m_hadrons) { particle.Rotate(rot_mat); }
    for(auto &particle : m_leptons) { particle.Rotate(rot_mat); }
}
