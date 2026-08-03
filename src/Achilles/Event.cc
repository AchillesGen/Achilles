#include "Achilles/Event.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::Event;

Event::Event(ProcessInfo& pi,std::shared_ptr<Nucleus> nuc=nullptr,double vwgt=0.0)
    : m_processInfo{pi}, m_nuc{nuc}, m_wgt{std::move(vwgt)} {
    m_hadrons = nuc->GenerateConfig();
}

Event &Event::operator=(const Event &other) {
    if(this == &other) return *this;
	m_processInfo = other.m_processInfo;
    m_nuc = other.m_nuc;
    m_remnant = other.m_remnant;
    m_wgt = other.m_wgt;
    leptonsIn = other.leptonsIn;
    leptonsOut = other.leptonsOut;
    hadronsIn = other.hadronsIn;
    hadronsOut = other.hadronsOut;
	spectators = other.spectators;
    m_history = other.m_history;
    flux = other.flux;
    m_process_id = other.m_process_id;
    return *this;
}

void Event::Finalize() {
    size_t nA = 0, nZ = 0;
    for(const Particle& part:allHadrons()) {
        if(part.IsExternal()) continue;
        if(part.ID() == PID::proton()) nZ++;
        nA++;
    }
    m_remnant = NuclearRemnant(nA, nZ);
}

void Event::Display() const {
    spdlog::trace("Leptons:");
    size_t idx = 0;
    for(const auto &particle: leptonsIn) { spdlog::trace("\t{}: {}", ++idx, particle); }
    for(const auto &particle: leptonsOut) { spdlog::trace("\t{}: {}", ++idx, particle); }
    spdlog::trace("Hadrons:");
    idx = 0;
    for(const auto &particle: hadronsIn) { spdlog::trace("\t{}: {}", ++idx, particle); }
    for(const auto &particle: hadronsOut) { spdlog::trace("\t{}: {}", ++idx, particle); }
    spdlog::trace("Weight: {}", Weight());
}

inline achilles::vParticles concatenate(std::vector<achilles::vParticles> lists) {
	achilles::vParticles result;
	for(achilles::vParticles list:lists)
		result.insert(result.end(),list.begin(),list.end());
	return result;
}

achilles::vParticles Event::Particles() const {
	return concatenate({leptonsIn,leptonsOut,hadronsIn,hadronsOut,spectators});
}

achilles::vParticles Event::allHadrons() const {
	return concatenate({hadronsIn,hadronsOut});
}

achilles::crefParticles Event::Protons(ParticleStatus status) const {
	vParticles hadrons=allHadrons();
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::proton(); };
        return FilterParticles(hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::proton() && p.Status() == status;
    };
    return FilterParticles(hadrons, func);
}

achilles::crefParticles Event::Neutrons(ParticleStatus status) const {
	vParticles hadrons=allHadrons();
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::neutron(); };
        return FilterParticles(hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::neutron() && p.Status() == status;
    };
    return FilterParticles(hadrons, func);
}

achilles::crefParticles Event::Pions(ParticleStatus status) const {
	vParticles hadrons=allHadrons();
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) {
            return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0());
        };
        return FilterParticles(hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0()) &&
               p.Status() == status;
    };
    return FilterParticles(hadrons, func);
}

achilles::refParticles Event::Protons(ParticleStatus status) {
	vParticles hadrons=allHadrons();
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::proton(); };
        return FilterParticles(hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::proton() && p.Status() == status;
    };
    return FilterParticles(hadrons, func);
}

achilles::refParticles Event::Neutrons(ParticleStatus status) {
	vParticles hadrons=allHadrons();
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) { return p.ID() == PID::neutron(); };
        return FilterParticles(hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return p.ID() == PID::neutron() && p.Status() == status;
    };
    return FilterParticles(hadrons, func);
}

achilles::refParticles Event::Pions(ParticleStatus status) {
	vParticles hadrons=allHadrons();
    if(status == ParticleStatus::any) {
        auto func = [](const Particle &p) {
            return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0());
        };
        return FilterParticles(hadrons, func);
    }
    auto func = [status](const Particle &p) {
        return (p.ID() == PID::pionp() || p.ID() == -PID::pionp() || p.ID() == PID::pion0()) &&
               p.Status() == status;
    };
    return FilterParticles(hadrons, func);
}

void Event::Rotate(const std::array<double, 9> &rot_mat) {
    for(Particle& particle: allHadrons()) { particle.Rotate(rot_mat); }
    for(Particle& particle: allHadrons()) { particle.Rotate(rot_mat); }
}
