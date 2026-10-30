// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVENT_HH
#define EVENT_HH

#include <algorithm>
#include <memory>
#include <vector>

#include "fmt/format.h"

#include "Achilles/Achilles.hh"
#include "Achilles/EventVertex.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Particle.hh"

namespace achilles {

class PID;
class FourVector;
class Particle;
class Beam;
class NuclearModel;
class HistoryVisitor;

using vParticles = std::vector<Particle>;
using refParticles = std::vector<std::reference_wrapper<Particle>>;
using crefParticles = std::vector<std::reference_wrapper<const Particle>>;


class Event {
  public:
    Event() = default;
	Event(ProcessInfo& pi,std::shared_ptr<Nucleus> nuc,double vwgt=0.0);
	Event(const Event&);
    Event &operator=(const Event &);
    MOCK ~Event() = default;

    void Finalize();

	/// Returns the particle representing the remnant that is currently in
	/// the event's final state (no end vertex), or nullptr if none exists.
	Particle* Remnant();
	/// Returns a pointer to the particle matching the argument
	/// within the event's final state, or nullptr if none exists.
	Particle* findFinalStateParticle(const Particle&);

    MOCK const std::shared_ptr<Nucleus> &CurrentNucleus() const { return m_nuc; }
    MOCK std::shared_ptr<Nucleus> &CurrentNucleus() { return m_nuc; }
	const ProcessInfo& processInfo() const { return m_processInfo; }

    const double &Flux() const { return flux; }
    double &Flux() { return flux; }
    MOCK const double &Weight() const { return m_wgt; }
    MOCK double &Weight() { return m_wgt; }
    void Rotate(const std::array<double, 9> &);
    void Display() const;

	vParticles& Particles() { return hist_particles; }

	/// Resets the event to its starting configuration
	/// so it can be passed through the mappers again.
	/// Intended for training/optimization purposes.
	void reset() {
		vertices.clear();
		hist_particles.clear();
	}

	refParticles getAllOfType(PID,ParticleStatus=ParticleStatus::any);
	crefParticles getAllOfType(PID,ParticleStatus=ParticleStatus::any) const;
	refParticles getAllOfType(refParticles&,PID,ParticleStatus=ParticleStatus::any);
	crefParticles getAllOfType(const crefParticles&,PID,ParticleStatus=ParticleStatus::any) const;

    bool operator==(const Event &other) const {
        return m_nuc==other.m_nuc
				&&vertices==other.vertices
				&&hist_particles==other.hist_particles;
    }

    int &ProcessId() { return m_process_id; }
    const int &ProcessId() const { return m_process_id; }

	/// Returns the new index of the created vertex
	size_t addVertex(EventVertex::StatusCode,const ThreeVector& = {});
	void addParticle(const Particle&,size_t,size_t);
	void addParticle(Particle&& part,size_t start_vertex,size_t end_vertex) {
		addParticle(part,start_vertex,end_vertex);
	}
	void addParticle(const Particle&,EventVertex*,EventVertex*);

	EventVertex* vertex(size_t) const;
	size_t nVertices() const { return vertices.size(); }
	/// Accepts a Vertex Index, returns other Vertex Indices
	std::vector<size_t> ChildVertices(size_t) const;
	/// Accepts a Vertex Index, returns other Vertex Indices
	std::vector<size_t> ParentVertices(size_t) const;
	EventVertex* uniqueVertex(EventVertex::StatusCode) const;
	/// Current implementation: primary is always at 0
	EventVertex* Primary() const { return vertices[0].get(); }

	/// Returns all particles that don't currently start on a vertex
	refParticles allInitialStateParticles();
	crefParticles allInitialStateParticles() const;
	/// Returns all particles that don't currently end on a vertex
	refParticles allFinalStateParticles();
	crefParticles allFinalStateParticles() const;

	/// Collects all current final-state hadrons into a list.
	/// Particles in this list are copies of the ones in the history.
	vParticles InteractableHadrons();

	/// A separate "slot" that is periodically updated to
	/// preserve compatibility w/ the various CrossSection classes.
	/// This is (hopefully) a temporary solution, and this field
	/// should not be treated as up-to-date in general.
	vParticles& CascadeHadrons() { return cascade_hadrons; }

	refParticles particlesIn(size_t);
	crefParticles particlesIn(size_t) const;
	refParticles particlesIn(size_t,PID,ParticleStatus=ParticleStatus::any);
	crefParticles particlesIn(size_t,PID,ParticleStatus=ParticleStatus::any) const;
	refParticles particlesOut(size_t);
	crefParticles particlesOut(size_t) const;
	refParticles particlesOut(size_t,PID,ParticleStatus=ParticleStatus::any);
	crefParticles particlesOut(size_t,PID,ParticleStatus=ParticleStatus::any) const;

	/// When the size_t argument is 0,
	/// these methods correspond to the ProcessInfo
	const Particle& leptonIn(size_t=0) const;
	crefParticles hadronsIn(size_t=0) const;
	crefParticles leptonsOut(size_t=0) const;
	crefParticles hadronsOut(size_t=0) const;
	crefParticles spectators() const;

	static vParticles MakeCopies(const crefParticles&);

	// Traversal
	void WalkHistory(HistoryVisitor&);

  private:
    // Helper functions
    refParticles FilterParticles(vParticles& particles,std::function<bool(const Particle&)> pred) {
        refParticles result;
		for(Particle& p:particles)
			if(pred(p))
				result.push_back(std::reference_wrapper(p));
        return result;
    }
    crefParticles FilterParticles(const vParticles& particles,std::function<bool(const Particle&)> pred) const {
        crefParticles result;
		for(const Particle& p:particles)
			if(pred(p))
				result.push_back(std::reference_wrapper(p));
        return result;
    }
    crefParticles FilterParticles(const crefParticles& particles,std::function<bool(const Particle&)> pred) const {
        crefParticles result;
		for(const Particle& p:particles)
			if(pred(p))
				result.push_back(std::reference_wrapper(p));
        return result;
    }
    refParticles FilterParticles(refParticles& particles,std::function<bool(const Particle&)> pred) {
        refParticles result;
		for(Particle& p:particles)
			if(pred(p))
				result.push_back(std::reference_wrapper(p));
        return result;
    }
    Particle MakeNuclearRemnant(size_t nA, size_t nZ, const FourVector &mom) const;

    // Variables
	ProcessInfo m_processInfo;
    std::shared_ptr<Nucleus> m_nuc;
    double m_wgt{};
	std::vector<std::unique_ptr<EventVertex>> vertices{};
    vParticles hist_particles{};
	vParticles cascade_hadrons{};
    double flux{};
    int m_process_id{};
};

class HistoryVisitor {
  public:
    virtual ~HistoryVisitor() = default;
    virtual void visit(Event&,size_t) = 0;
};

struct PrintVisitor : HistoryVisitor {
	std::string data;
	void visit(Event& event,size_t vertexIdx) {
		refParticles in=event.particlesIn(vertexIdx),out=event.particlesOut(vertexIdx);
		data+=fmt::format("Vertex({}, {{{}}} -> {{{}}})\n", achilles::ToString(event.vertex(vertexIdx)->Status()),
				fmt::join(in.begin(), in.end(), ", "),fmt::join(out.begin(), out.end(), ", "));
	}
};

struct PIDLocator : HistoryVisitor {
    std::vector<Particle> particles;
    PID pid;
    int direction;
    PIDLocator(PID _pid, int dir) : pid(_pid), direction(dir) {}
    void visit(Event& event,size_t vertexIdx) {
        if(direction >= 0) {
            for(Particle& particle:event.particlesOut(vertexIdx)) {
                if(particle.ID() == pid) particles.push_back(particle);
            }
        } else {
            for(Particle& particle:event.particlesIn(vertexIdx)) {
                if(particle.ID() == pid) particles.push_back(particle);
            }
        }
    }
};

} // namespace achilles

template <> struct fmt::formatter<achilles::Event> {
    char presentation = 'e';
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        // Parse the presentation format and store it in the formatter:
        auto it = ctx.begin(), end = ctx.end();
        if(it != end && (*it == 'f' || *it == 'e')) presentation = *it++;

        // Check if reached the end of the range:
        if(it != end && *it != '}') throw format_error("Invalid format");

        // Return an iterator past the end of the parsed range:
        return it;
    }

    auto format(const achilles::Event&, format_context &ctx) const
        -> format_context::iterator {
        // ctx.out() is an output iterator to write to
        return format_to(ctx.out(),"Event. (TODO: Implement Formatted Output)");
    }
};

#endif