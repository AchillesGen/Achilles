#include "Achilles/Event.hh"
#include "Achilles/Beams.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::Event;

Event::Event(std::shared_ptr<Nucleus> nuc, std::vector<FourVector> mom, double vwgt)
    : m_nuc{std::move(nuc)}, m_mom{std::move(mom)}, m_wgt{std::move(vwgt)} {
    m_nuc->GenerateConfig();
    m_me.resize(m_nuc->NNucleons());
}

void Event::Finalize() {
    size_t nA = 0, nZ = 0;
    for(auto it = m_nuc->Nucleons().begin(); it != m_nuc->Nucleons().end();) {
        if(it->Status() == ParticleStatus::background ||
           it->Status() == ParticleStatus::propagating) {
            if(it->ID() == PID::proton()) nZ++;
            nA++;
            it = m_nuc->Nucleons().erase(it);
        } else {
            ++it;
        }
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
    result.insert(result.end(), Hadrons().begin(), Hadrons().end());
    result.insert(result.end(), m_leptons.begin(), m_leptons.end());
    return result;
}

const achilles::vParticles &Event::Hadrons() const {
    return m_nuc->Nucleons();
}

achilles::vParticles &Event::Hadrons() {
    return m_nuc->Nucleons();
}

void Event::Rotate(const std::array<double, 9> &rot_mat) {
    for(auto &particle : m_nuc->Nucleons()) { particle.Rotate(rot_mat); }
    for(auto &particle : m_leptons) { particle.Rotate(rot_mat); }
}
