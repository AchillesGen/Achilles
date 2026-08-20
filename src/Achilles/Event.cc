// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Event.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::Event;

Event::Event(std::shared_ptr<Nucleus> nuc, std::vector<FourVector> mom, double vwgt)
    : m_nuc{nuc}, m_mom{std::move(mom)}, m_wgt{std::move(vwgt)} {
    m_hadrons = nuc->GenerateConfig();
}

Event::Event(const Event &other) {
    m_nuc = other.m_nuc;
    m_remnant = other.m_remnant;
    m_mom = other.m_mom;
    m_wgt = other.m_wgt;
    m_leptons = other.m_leptons;
    m_hadrons = other.m_hadrons;
    m_history = other.m_history;
    flux = other.flux;
    m_process_id = other.m_process_id;
}

Event &Event::operator=(const Event &other) {
    if(this == &other) return *this;
    m_nuc = other.m_nuc;
    m_remnant = other.m_remnant;
    m_mom = other.m_mom;
    m_wgt = other.m_wgt;
    m_leptons = other.m_leptons;
    m_hadrons = other.m_hadrons;
    m_history = other.m_history;
    flux = other.flux;
    m_process_id = other.m_process_id;
    return *this;
}

void Event::Finalize() {
    spdlog::debug("Finalizing the event");
    if(m_remnant.ID() == PID::undefined()) return;

    // Every residue is created as intermediate; the one still standing at the end is the
    // actual remnant. Promoting it here is what lets a reader pick it out of the chain by
    // status instead of walking the graph for an unconsumed residue.
    m_remnant.Status() = ParticleStatus::residue;
    m_history.UpdateStatuses({m_remnant});

    // The remnant is built up through the history rather than derived here. Cross-check it
    // against global four-momentum balance, which must agree if every vertex conserved.
    if(spdlog::get_level() <= spdlog::level::debug && m_nuc && !m_leptons.empty()) {
        FourVector mom = m_leptons[0].Momentum() + m_nuc->InitParticle().Momentum();
        for(size_t i = 1; i < m_leptons.size(); ++i) { mom -= m_leptons[i].Momentum(); }
        for(const auto &part : m_hadrons) {
            if(part.IsFinal()) mom -= part.Momentum();
        }
        const auto diff = mom - m_remnant.Momentum();
        if(std::abs(diff.E()) > 1e-6 || diff.Vec3().Magnitude() > 1e-6)
            spdlog::debug("Event: remnant disagrees with momentum balance by {}", diff);
    }
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

achilles::crefParticles Event::Pions(ParticleStatus status) const {
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) {
            return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0());
        };
        return FilterParticles(m_hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0()) &&
               p.Status() == status;
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

achilles::refParticles Event::Pions(ParticleStatus status) {
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) {
            return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0());
        };
        return FilterParticles(m_hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0()) &&
               p.Status() == status;
    };
    return FilterParticles(m_hadrons, func);
}

void Event::Rotate(const std::array<double, 9> &rot_mat) {
    for(auto &particle : m_hadrons) { particle.Rotate(rot_mat); }
    for(auto &particle : m_leptons) { particle.Rotate(rot_mat); }
}

achilles::Particle Event::MakeNuclearRemnant(size_t nA, size_t nZ, const FourVector &mom) const {
    PID pid = PID::MakeNucleus(static_cast<int>(nZ), static_cast<int>(nA));
    return Particle(pid, mom, {}, ParticleStatus::residue);
}
