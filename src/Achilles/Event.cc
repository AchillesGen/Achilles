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
    for(const Particle& particle: leptonsIn) { spdlog::trace("\t{}: {}", ++idx, particle); }
    for(const Particle& particle: leptonsOut) { spdlog::trace("\t{}: {}", ++idx, particle); }
    spdlog::trace("Mapper Hadrons:");
    idx = 0;
    for(const Particle& particle: hadronsIn) { spdlog::trace("\t{}: {}", ++idx, particle); }
    for(const Particle& particle: hadronsOut) { spdlog::trace("\t{}: {}", ++idx, particle); }
	spdlog::trace("Cascade Hadrons:");
    idx = 0;
    for(const Particle& particle: nucleus_hadrons) { spdlog::trace("\t{}: {}", ++idx, particle); }
    spdlog::trace("Weight: {}", Weight());
}

achilles::ptrParticles Event::getAllOfType(vParticles& list,PID pid,ParticleStatus status) {
    if(status == ParticleStatus::any) {
		auto func = [pid](const Particle& p) { return p.ID()==pid; };
		return FilterPointers(list,func);
	}
	auto func = [pid,status](const Particle& p) { return p.ID()==pid&&p.Status()==status; };
	return FilterPointers(list,func);
}

achilles::ptrParticles Event::allParticles() {
	ptrParticles result;
	for(Particle& p:leptonsIn)
		result.push_back(&p);
	for(Particle& p:leptonsOut)
		result.push_back(&p);
	if(hadrons_setup) {
		for(Particle& p:nucleus_hadrons)
			result.push_back(&p);
	} else {
		for(Particle& p:hadronsIn)
			result.push_back(&p);
		for(Particle& p:hadronsOut)
			result.push_back(&p);
		for(Particle& p:spectators)
			result.push_back(&p);
	}
	return result;
}
achilles::cptrParticles Event::allParticles() const {
	cptrParticles result;
	for(const Particle& p:leptonsIn)
		result.push_back(&p);
	for(const Particle& p:leptonsOut)
		result.push_back(&p);
	if(hadrons_setup) {
		for(const Particle& p:nucleus_hadrons)
			result.push_back(&p);
	} else {
		for(const Particle& p:hadronsIn)
			result.push_back(&p);
		for(const Particle& p:hadronsOut)
			result.push_back(&p);
		for(const Particle& p:spectators)
			result.push_back(&p);
	}
	return result;
}

void Event::Rotate(const std::array<double, 9> &rot_mat) {
    for(Particle* particle:allParticles()) { particle->Rotate(rot_mat); }
}

void Event::assignParticleDetails(vParticles& particleSource,PID pid) {
	if(particleSource.empty())
		return;
	ptrParticles sources=getAllOfType(particleSource,pid);
	if(sources.empty())
		return;
	ptrParticles candidates=getAllOfType(nucleus_hadrons,pid,ParticleStatus::background);
	std::vector<size_t> targets=Random::Instance().SampleIndices(candidates.size(),sources.size());
	for(size_t i=0;i<sources.size();i++) {
		// Target gets Status and Momentum of mapped initial state
		candidates[targets[i]]->Status()=sources[i]->Status();
		candidates[targets[i]]->Momentum()=sources[i]->Momentum();
		// Initial state gets Position of target, for finalstate calc later
		sources[i]->Position()=candidates[targets[i]]->Position();
	}
}

void Event::SetupHadrons() {
	// Coherent Scattering case
	if(ParticleInfo(m_processInfo->m_hadronic.first[0]).IsNucleus())
		return;

	assignParticleDetails(hadronsIn,PID::proton());
	assignParticleDetails(hadronsIn,PID::neutron());
	assignParticleDetails(spectators,PID::proton());
	assignParticleDetails(spectators,PID::neutron());

    // Initialize final state hadrons
    // TODO: Handle propagating deltas
    // TODO: Handle selecting position for things like MEC+pion production
    size_t cur_idx = 0;
    ThreeVector position;
    for(size_t i = 0; i < hadronsOut.size(); ++i) {
        if(ParticleInfo(m_processInfo->m_hadronic.second[i]).IsBaryon())
            position = hadronsIn[cur_idx++].Position();
        Particle& part=hadronsOut[i];
        part.Status() = ParticleStatus::propagating;
        part.Position() = position;
		nucleus_hadrons.push_back(part);
    }
	hadrons_setup=true;
}