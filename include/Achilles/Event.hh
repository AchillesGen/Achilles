#ifndef EVENT_HH
#define EVENT_HH

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "fmt/format.h"

#include "Achilles/Achilles.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/NuclearRemnant.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Particle.hh"

namespace achilles {

class PID;
class FourVector;
class Particle;
class Beam;
class NuclearModel;

using vParticles = std::vector<Particle>;
using ptrParticles = std::vector<Particle*>;
using cptrParticles = std::vector<const Particle*>;
using vMomentum = std::vector<FourVector>;
using refParticles = std::vector<std::reference_wrapper<Particle>>;
using crefParticles = std::vector<std::reference_wrapper<const Particle>>;

class Event {
  public:
    Event() = default;
	Event(const Event&);
	Event(ProcessInfo& pi,std::shared_ptr<Nucleus> nuc,double vwgt=0.0):
			m_processInfo{std::make_shared<ProcessInfo>(pi)}, m_nuc{nuc}, m_wgt{std::move(vwgt)} {
				spdlog::trace("Event Process: {}",pi);
			}
    Event &operator=(const Event &);
    MOCK ~Event() = default;

    void Finalize();

    MOCK const NuclearRemnant &Remnant() const { return m_remnant; }

    MOCK const std::shared_ptr<Nucleus> &CurrentNucleus() const { return m_nuc; }
    MOCK std::shared_ptr<Nucleus> &CurrentNucleus() { return m_nuc; }
	ProcessInfo& processInfo() { return *m_processInfo; }

    const double &Flux() const { return flux; }
    double &Flux() { return flux; }
    MOCK const double &Weight() const { return m_wgt; }
    MOCK double &Weight() { return m_wgt; }
    void Rotate(const std::array<double, 9> &);
    void Display() const;

    MOCK const vParticles &NucleusHadrons() const { return nucleus_hadrons; }
    MOCK vParticles &NucleusHadrons() { return nucleus_hadrons; }
    MOCK const vParticles &LeptonsIn() const { return leptonsIn; }
    MOCK vParticles &LeptonsIn() { return leptonsIn; }
    MOCK const vParticles &LeptonsOut() const { return leptonsOut; }
    MOCK vParticles &LeptonsOut() { return leptonsOut; }
    MOCK const vParticles &HadronsIn() const { return hadronsIn; }
    MOCK vParticles &HadronsIn() { return hadronsIn; }
    MOCK const vParticles &HadronsOut() const { return hadronsOut; }
    MOCK vParticles &HadronsOut() { return hadronsOut; }
    MOCK const vParticles &Spectators() const { return spectators; }
    MOCK vParticles &Spectators() { return spectators; }

	/// Resets the event to its starting configuration
	/// so it can be passed through the mappers again.
	/// Intended for training/optimization purposes.
	void reset() {
		leptonsIn.clear();
		hadronsIn.clear();
		leptonsOut.clear();
		hadronsOut.clear();
		spectators.clear();
		hadrons_setup=false;
	}
	void addLeptonIn(FourVector momentum,ParticleStatus status=ParticleStatus::initial_state) {
		if(leptonsIn.size()>=1)
			throw std::runtime_error("Event::addLeptonIn(): Particle count exceeds process specifications");
		spdlog::trace("Creating Lepton-In ({}, {})",m_processInfo->m_leptonic.first,status);
		leptonsIn.push_back(Particle(m_processInfo->m_leptonic.first,momentum,{},status));
	}
	void addLeptonOut(FourVector momentum,ParticleStatus status=ParticleStatus::final_state) {
		if(leptonsOut.size()>=m_processInfo->m_leptonic.second.size())
			throw std::runtime_error("Event::addLeptonOut(): Particle count exceeds process specifications");
		spdlog::trace("Creating Lepton-Out ({}, {})",m_processInfo->m_leptonic.second[leptonsOut.size()],status);
		leptonsOut.push_back(Particle(m_processInfo->m_leptonic.second[leptonsOut.size()],momentum,{},status));
	}
	void addHadronIn(FourVector momentum,ParticleStatus status=ParticleStatus::initial_state) {
		if(hadronsIn.size()>=m_processInfo->m_hadronic.first.size())
			throw std::runtime_error("Event::addHadronIn(): Particle count exceeds process specifications");
		spdlog::trace("Creating Hadron-In ({}, {})",m_processInfo->m_hadronic.first[hadronsIn.size()],status);
		hadronsIn.push_back(Particle(m_processInfo->m_hadronic.first[hadronsIn.size()],momentum,{},status));
	}
	void addHadronOut(FourVector momentum,ParticleStatus status=ParticleStatus::final_state) {
		if(hadronsOut.size()>=m_processInfo->m_hadronic.second.size())
			throw std::runtime_error("Event::addHadronOut(): Particle count exceeds process specifications");
		spdlog::trace("Creating Hadron-Out ({} ,{})",m_processInfo->m_hadronic.second[hadronsOut.size()],status);
		hadronsOut.push_back(Particle(m_processInfo->m_hadronic.second[hadronsOut.size()],momentum,{},status));
	}
	void addSpectator(FourVector momentum,ParticleStatus status=ParticleStatus::spectator) {
		if(spectators.size()>=m_processInfo->m_spectator.size())
			throw std::runtime_error("Event::addSpectator(): Particle count exceeds process specifications");
		spdlog::trace("Creating Spectator ({}, {})",m_processInfo->m_spectator[spectators.size()],status);
		spectators.push_back(Particle(m_processInfo->m_spectator[spectators.size()],momentum,{},status));
	}
	void addAutoOutgoing(FourVector momentum,ParticleStatus status=ParticleStatus::final_state) {
		if(leptonsOut.size()<m_processInfo->m_leptonic.second.size())
			addLeptonOut(momentum,status);
		else
			addHadronOut(momentum,status);
	}

	void SetupHadrons();

	ptrParticles getAllOfType(vParticles&,PID,ParticleStatus=ParticleStatus::any);
	ptrParticles allParticles();
	cptrParticles allParticles() const;

    MOCK const EventHistory &History() const { return m_history; }
    EventHistory &History() { return m_history; }

    bool operator==(const Event &other) const {
        return m_nuc == other.m_nuc && m_remnant == other.m_remnant
				&&nucleus_hadrons==other.nucleus_hadrons
				&&leptonsIn==other.leptonsIn && leptonsOut==other.leptonsOut
				&&hadronsIn==other.hadronsIn && hadronsOut==other.hadronsOut;
    }

    int &ProcessId() { return m_process_id; }
    const int &ProcessId() const { return m_process_id; }

  private:
    // Helper functions
    template <class UnaryPred> ptrParticles FilterPointers(vParticles& particles, UnaryPred pred) {
        ptrParticles result;
		for(Particle& p:particles)
			if(pred(p))
				result.push_back(&p);
        return result;
    }
    template <class UnaryPred>
    crefParticles FilterParticles(const crefParticles& particles, UnaryPred pred) const {
        crefParticles result;
        std::copy_if(particles.begin(), particles.end(), std::back_inserter(result), pred);
        return result;
    }
    template <class UnaryPred> refParticles FilterParticles(refParticles& particles, UnaryPred pred) {
        refParticles result;
        std::copy_if(particles.begin(), particles.end(), std::back_inserter(result), pred);
        return result;
    }

	vParticles concatenate(std::vector<vParticles> lists) const {
		vParticles result;
		for(vParticles list:lists)
			std::copy(list.begin(), list.end(), std::back_inserter(result));
		return result;
	}
	/// Takes all particles of the given PID from the given list of particles,
	/// and assigns their Status and Momentum to randomly-selected particles of
	/// the same type in the given event's "nucleus_hadrons" list.
    void assignParticleDetails(vParticles&,PID);

    // Variables
	std::shared_ptr<ProcessInfo> m_processInfo;
    std::shared_ptr<Nucleus> m_nuc;
    NuclearRemnant m_remnant{};
    //vMomentum m_mom{};
    double m_wgt{};
    vParticles nucleus_hadrons{}, leptonsIn{}, leptonsOut{}, hadronsIn{}, hadronsOut{}, spectators{};
    EventHistory m_history{};
    double flux{};
    int m_process_id{};
	bool hadrons_setup=false;
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