#include "Achilles/Event.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

using achilles::Event;

Event::Event(const Event& other) {
	m_processInfo = other.m_processInfo;
    m_nuc = other.m_nuc;
    m_remnant = other.m_remnant;
    m_wgt = other.m_wgt;
	nucleus_hadrons = other.nucleus_hadrons;
    leptonsIn = other.leptonsIn;
    leptonsOut = other.leptonsOut;
    hadronsIn = other.hadronsIn;
    hadronsOut = other.hadronsOut;
	spectators = other.spectators;
    m_history = other.m_history;
    flux = other.flux;
    m_process_id = other.m_process_id;
}

Event &Event::operator=(const Event &other) {
    if(this == &other) return *this;
	m_processInfo = other.m_processInfo;
    m_nuc = other.m_nuc;
    m_remnant = other.m_remnant;
    m_wgt = other.m_wgt;
	nucleus_hadrons = other.nucleus_hadrons;
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
    for(const Particle& part:NucleusHadrons()) {
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

achilles::refParticles Event::getAllOfType(vParticles list,PID pid,ParticleStatus status) {
    if(status == ParticleStatus::any) {
		auto func = [pid](const Particle& p) { return p.ID()==pid; };
		return FilterParticles(list,func);
	}
	auto func = [pid,status](const Particle& p) { return p.ID()==pid&&p.Status()==status; };
	return FilterParticles(list,func);
}
achilles::refParticles Event::getAllOfType(refParticles list,PID pid,ParticleStatus status) {
    if(status == ParticleStatus::any) {
		auto func = [pid](const Particle& p) { return p.ID()==pid; };
		return FilterParticles(list,func);
	}
	auto func = [pid,status](const Particle& p) { return p.ID()==pid&&p.Status()==status; };
	return FilterParticles(list,func);
}

achilles::crefParticles Event::allParticles() const {
	return concatenate({nucleus_hadrons,leptonsIn,leptonsOut,hadronsIn,hadronsOut,spectators});
}
achilles::refParticles Event::allParticles() {
	return concatenate({nucleus_hadrons,leptonsIn,leptonsOut,hadronsIn,hadronsOut,spectators});
}

achilles::crefParticles Event::allHadrons() const {
	return concatenate({nucleus_hadrons,hadronsIn,hadronsOut});
}
achilles::refParticles Event::allHadrons() {
	return concatenate({nucleus_hadrons,hadronsIn,hadronsOut});
}

void Event::Rotate(const std::array<double, 9> &rot_mat) {
    for(Particle& particle:allParticles()) { particle.Rotate(rot_mat); }
}