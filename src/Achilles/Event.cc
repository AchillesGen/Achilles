// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Event.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

#include <queue>

using achilles::Event;
using achilles::EventVertex;
using achilles::Particle;

Event::Event(ProcessInfo& pi,std::shared_ptr<Nucleus> nuc,double vwgt):
			m_processInfo{pi}, m_nuc{nuc}, m_wgt{std::move(vwgt)} {
	// Primary vertex is always at position 0
	addVertex(EventVertex::StatusCode::primary);
	spdlog::trace("Event Process: {}",pi);
}

Event::Event(const Event& other) {
	m_processInfo = other.m_processInfo;
    m_nuc = other.m_nuc;
    m_wgt = other.m_wgt;
	hist_particles = other.hist_particles;
	cascade_hadrons = other.cascade_hadrons;
    flux = other.flux;
    m_process_id = other.m_process_id;

    vertices.reserve(other.vertices.size());
    for(const auto &elm : other.vertices) {
        vertices.push_back(std::make_unique<EventVertex>(*elm));
    }
}

Event &Event::operator=(const Event &other) {
    if(this == &other) return *this;
	m_processInfo = other.m_processInfo;
    m_nuc = other.m_nuc;
    m_wgt = other.m_wgt;
	hist_particles = other.hist_particles;
	cascade_hadrons = other.cascade_hadrons;
    flux = other.flux;
    m_process_id = other.m_process_id;

    vertices.clear();
    vertices.reserve(other.vertices.size());
    for(const auto &elm : other.vertices) {
        vertices.push_back(std::make_unique<EventVertex>(*elm));
    }

    return *this;
}

void Event::Finalize() {
    // Every residue is created as intermediate; the one still standing at the end is the
    // actual remnant. Promoting it here is what lets a reader pick it out of the chain by
    // status instead of walking the graph for an unconsumed residue.
	if(Particle* rem=Remnant())
		rem->Status()=ParticleStatus::residue;

    // The remnant is built up through the history rather than derived here. Cross-check it
    // against global four-momentum balance, which must agree if every vertex conserved.

	// This code was equivalent to just checking global momentum conservation,
	// so I made that more explicit -Hayden
	if(spdlog::get_level() <= spdlog::level::debug) {
		spdlog::debug("Ensuring Four-Momentum Conservation...");
		FourVector diff;
		for(Particle& part:hist_particles) {
			if(part.isInInitialState())
				diff-=part.Momentum();
			if(part.isInFinalState())
				diff+=part.Momentum();
		}
		if(std::abs(diff.E()) > 1e-6 || diff.Vec3().Magnitude() > 1e-6)
			spdlog::debug("Event: Total Four-Momentum not conserved, difference of {}", diff);
	}
}

achilles::NuclearRemnant* Event::Remnant() {
	for(Particle& part:hist_particles)
		if(part.isRemnant()&&part.isInFinalState())
			return static_cast<NuclearRemnant*>(&part);
	return nullptr;
}
Particle* Event::findFinalStateParticle(const Particle& match) {
	for(Particle& part:hist_particles)
		if(part.isInFinalState()&&part.equalsIgnoreVertices(match))
			return &part;
	return nullptr;
}

void Event::Display() const {
    spdlog::trace("Particles in History:");
    size_t idx = 0;
    for(const Particle& particle:hist_particles) { spdlog::trace("\t{}: {}", ++idx, particle); }
    spdlog::trace("Weight: {}", Weight());
}

achilles::refParticles Event::getAllOfType(PID pid,ParticleStatus status) {
    if(status == ParticleStatus::any) {
		auto func = [pid](const Particle& p) { return p.ID()==pid; };
		return FilterParticles(hist_particles,func);
	}
	auto func = [pid,status](const Particle& p) { return p.ID()==pid&&p.Status()==status; };
	return FilterParticles(hist_particles,func);
}
achilles::crefParticles Event::getAllOfType(PID pid,ParticleStatus status) const {
    if(status == ParticleStatus::any) {
		auto func = [pid](const Particle& p) { return p.ID()==pid; };
		return FilterParticles(hist_particles,func);
	}
	auto func = [pid,status](const Particle& p) { return p.ID()==pid&&p.Status()==status; };
	return FilterParticles(hist_particles,func);
}
achilles::refParticles Event::getAllOfType(refParticles& list,PID pid,ParticleStatus status) {
    if(status == ParticleStatus::any) {
		auto func = [pid](const Particle& p) { return p.ID()==pid; };
		return FilterParticles(list,func);
	}
	auto func = [pid,status](const Particle& p) { return p.ID()==pid&&p.Status()==status; };
	return FilterParticles(list,func);
}
achilles::crefParticles Event::getAllOfType(const crefParticles& list,PID pid,ParticleStatus status) const {
    if(status == ParticleStatus::any) {
		auto func = [pid](const Particle& p) { return p.ID()==pid; };
		return FilterParticles(list,func);
	}
	auto func = [pid,status](const Particle& p) { return p.ID()==pid&&p.Status()==status; };
	return FilterParticles(list,func);
}

void Event::Rotate(const std::array<double, 9> &rot_mat) {
    for(Particle& particle:hist_particles) { particle.Rotate(rot_mat); }
}

/*void Event::assignParticleDetails(vParticles& particleSource,PID pid) {
	if(particleSource.empty())
		return;
	refParticles sources=getAllOfType(particleSource,pid);
	if(sources.empty())
		return;
	refParticles candidates=getAllOfType(p_backend,pid,ParticleStatus::background);
	std::vector<size_t> targets=Random::Instance().SampleIndices(candidates.size(),sources.size());
	for(size_t i=0;i<sources.size();i++) {
		// Target gets Status and Momentum of mapped initial state
		candidates[targets[i]].get().Status()=sources[i].get().Status();
		candidates[targets[i]].get().Momentum()=sources[i].get().Momentum();
		// Initial state gets Position of target, for finalstate calc later
		sources[i].get().Position()=candidates[targets[i]].get().Position();
	}
}

void Event::SetupHadrons() {
	// Coherent Scattering case
	if(ParticleInfo(m_processInfo.m_hadronic.first[0]).IsNucleus())
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
}*/

size_t Event::addVertex(EventVertex::StatusCode status,const ThreeVector& position) {
	size_t idx=vertices.size();
    vertices.push_back(std::make_unique<EventVertex>(idx,status,position));
	return idx;
}

void Event::addParticle(const Particle& part,size_t start_vertex,size_t end_vertex) {
	hist_particles.push_back(part);
	hist_particles.back().startVertex()=start_vertex;
	hist_particles.back().endVertex()=end_vertex;
}
void Event::addParticle(const Particle& part,EventVertex* start,EventVertex* end) {
	addParticle(part,start?start->Index():Particle::OPEN_END,end?end->Index():Particle::OPEN_END);
}

std::vector<size_t> Event::ChildVertices(size_t vertexIdx) const {
	std::vector<size_t> result;
	for(const Particle& part:hist_particles) {
		// If the particle starts on the given vertex and ends on another valid vertex,
		// the end is included as a child of the start
		if(part.startVertex()==vertexIdx&&part.endVertex()!=Particle::OPEN_END)
			result.push_back(part.endVertex());
	}
	return result;
}

std::vector<size_t> Event::ParentVertices(size_t vertexIdx) const {
	std::vector<size_t> result;
	for(const Particle& part:hist_particles) {
		// If the particle ends on the given vertex and starts on another valid vertex,
		// the start is included as a parent of the end
		if(part.endVertex()==vertexIdx&&part.startVertex()!=Particle::OPEN_END)
			result.push_back(part.startVertex());
	}
	return result;
}

EventVertex* Event::vertex(size_t idx) const {
	if(idx<vertices.size())
		return vertices[idx].get();
	return nullptr;
}

EventVertex* Event::uniqueVertex(EventVertex::StatusCode status) const {
	EventVertex* result=nullptr;
	for(const std::unique_ptr<EventVertex>& node:vertices)
		if(node->Status()==status) {
			if(result==nullptr)
				result = node.get();
			else
				throw std::runtime_error(fmt::format("EventHistory: Only one {} vertex is allowed!",ToString(status)));
		}
	return result;
}

achilles::refParticles Event::allInitialStateParticles() {
	// Initial State particles list OPEN_END as their start vertex
	return particlesOut(Particle::OPEN_END);
}
achilles::crefParticles Event::allInitialStateParticles() const {
	return particlesOut(Particle::OPEN_END);
}

achilles::refParticles Event::allFinalStateParticles() {
	// Final State particles list OPEN_END as their end vertex
	return particlesIn(Particle::OPEN_END);
}
achilles::crefParticles Event::allFinalStateParticles() const {
	return particlesIn(Particle::OPEN_END);
}

achilles::vParticles Event::InteractableHadrons() {
	vParticles result;
	for(Particle& part:hist_particles) if(part.isInFinalState()) {
		if(part.isRemnant()) {
			for(Particle& sub:static_cast<NuclearRemnant*>(&part)->sub_parts)
				if(ParticleInfo(sub.ID()).IsHadron())
					result.push_back(part);
		} else if(ParticleInfo(part.ID()).IsHadron())
			result.push_back(part);
	}
	return result;
}

achilles::refParticles Event::particlesIn(size_t vertex) {
	refParticles result;
	for(Particle& part:hist_particles)
		if(part.endVertex()==vertex)
			result.push_back(std::reference_wrapper<Particle>(part));
	return result;
}
achilles::crefParticles Event::particlesIn(size_t vertex) const {
	crefParticles result;
	for(const Particle& part:hist_particles)
		if(part.endVertex()==vertex)
			result.push_back(std::reference_wrapper<const Particle>(part));
	return result;
}
achilles::refParticles Event::particlesIn(size_t vertex,PID pid,ParticleStatus status) {
	refParticles list=particlesIn(vertex);
	return getAllOfType(list,pid,status);
}
achilles::crefParticles Event::particlesIn(size_t vertex,PID pid,ParticleStatus status) const {
	crefParticles list=particlesIn(vertex);
	return getAllOfType(list,pid,status);
}
achilles::refParticles Event::particlesOut(size_t vertex) {
	refParticles result;
	for(Particle& part:hist_particles)
		if(part.startVertex()==vertex)
			result.push_back(std::reference_wrapper<Particle>(part));
	return result;
}
achilles::crefParticles Event::particlesOut(size_t vertex) const {
	crefParticles result;
	for(const Particle& part:hist_particles)
		if(part.startVertex()==vertex)
			result.push_back(std::reference_wrapper<const Particle>(part));
	return result;
}
achilles::refParticles Event::particlesOut(size_t vertex,PID pid,ParticleStatus status) {
	refParticles list=particlesOut(vertex);
	return getAllOfType(list,pid,status);
}
achilles::crefParticles Event::particlesOut(size_t vertex,PID pid,ParticleStatus status) const {
	crefParticles list=particlesOut(vertex);
	return getAllOfType(list,pid,status);
}

const Particle& Event::leptonIn(size_t vertex) const {
	return particlesIn(vertex,m_processInfo.m_leptonic.first)[0].get();
}
achilles::crefParticles Event::hadronsIn(size_t vertex) const {
	auto pred=[*this](const Particle& p) {
		for(PID pid:m_processInfo.m_hadronic.first)
			if(p.ID()==pid)
				return true;
		return false;
	};
	return FilterParticles(particlesIn(vertex),pred);
}
achilles::crefParticles Event::leptonsOut(size_t vertex) const {
	auto pred=[*this](const Particle& p) {
		for(PID pid:m_processInfo.m_leptonic.second)
			if(p.ID()==pid)
				return true;
		return false;
	};
	return FilterParticles(particlesOut(vertex),pred);
}
achilles::crefParticles Event::hadronsOut(size_t vertex) const {
	auto pred=[*this](const Particle& p) {
		for(PID pid:m_processInfo.m_hadronic.second)
			if(p.ID()==pid)
				return true;
		return false;
	};
	return FilterParticles(particlesOut(vertex),pred);
}
achilles::crefParticles Event::spectators() const {
	auto pred=[*this](const Particle& p) {
		return p.Status()==ParticleStatus::spectator;
	};
	return FilterParticles(hist_particles,pred);
}

achilles::vParticles Event::MakeCopies(const crefParticles& refs) {
	vParticles result;
	for(const Particle& p:refs)
		result.push_back(p);
	return result;
}

void Event::WalkHistory(achilles::HistoryVisitor &visitor) {
    std::vector<size_t> visited;
    std::queue<size_t> to_visit;

    // Start at Primary Interaction
	size_t current=Primary()->Index();
    to_visit.push(current);

    while(!to_visit.empty()) {
        current = to_visit.front();
        to_visit.pop();

        // Ensure it hasn't been visited
        if(std::find(visited.begin(), visited.end(), current) != visited.end()) continue;

        visitor.visit(*this,current);

        for(size_t& child:ChildVertices(current)) to_visit.push(child);

        for(size_t& parent:ParentVertices(current)) to_visit.push(parent);

        visited.push_back(current);
    }
}

achilles::Particle Event::MakeNuclearRemnant(size_t nA, size_t nZ, const FourVector &mom) const {
    PID pid = PID::MakeNucleus(static_cast<int>(nZ), static_cast<int>(nA));
    return Particle(pid, mom, {}, ParticleStatus::residue);
}