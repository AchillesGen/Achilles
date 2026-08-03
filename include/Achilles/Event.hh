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
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Particle.hh"

namespace achilles {

class PID;
class FourVector;
class Particle;
class Nucleus;
class Beam;
class NuclearModel;

using vParticles = std::vector<Particle>;
using vMomentum = std::vector<FourVector>;
using refParticles = std::vector<std::reference_wrapper<Particle>>;
using crefParticles = std::vector<std::reference_wrapper<const Particle>>;

class Event {
  public:
    Event() = default;
    Event(ProcessInfo&,std::shared_ptr<Nucleus> =nullptr, double=0.0);
    Event &operator=(const Event &);
    MOCK ~Event() = default;

    void Finalize();

    MOCK const NuclearRemnant &Remnant() const { return m_remnant; }

	//MOCK const vMomentum &Momentum() const { return m_mom; }
	//MOCK vMomentum &Momentum() { return m_mom; }

    MOCK const std::shared_ptr<Nucleus> &CurrentNucleus() const { return m_nuc; }
    MOCK std::shared_ptr<Nucleus> &CurrentNucleus() { return m_nuc; }
	void setNucleus(std::shared_ptr<Nucleus> nucleus) { m_nuc=nucleus; }
	ProcessInfo& processInfo() { return m_processInfo; }

    const double &Flux() const { return flux; }
    double &Flux() { return flux; }

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
    MOCK const double &Weight() const { return m_wgt; }
    MOCK double &Weight() { return m_wgt; }
    void Rotate(const std::array<double, 9> &);
    void Display() const;

	void addLeptonIn(FourVector momentum,ParticleStatus status=ParticleStatus::beam) {
		leptonsIn.push_back(Particle(m_processInfo.m_leptonic.first,momentum,{},status));
	}
	void addLeptonOut(FourVector momentum,ParticleStatus status=ParticleStatus::final_state) {
		leptonsOut.push_back(Particle(m_processInfo.m_leptonic.second[leptonsOut.size()],momentum,{},status));
	}
	void addHadronIn(FourVector momentum,ParticleStatus status=ParticleStatus::decayed) {
		hadronsIn.push_back(Particle(m_processInfo.m_hadronic.first[hadronsIn.size()],momentum,{},status));
	}
	void addHadronOut(FourVector momentum,ParticleStatus status=ParticleStatus::final_state) {
		hadronsOut.push_back(Particle(m_processInfo.m_hadronic.second[hadronsOut.size()],momentum,{},status));
	}
	void addSpectator(FourVector momentum,ParticleStatus status=ParticleStatus::spectator) {
		spectators.push_back(Particle(m_processInfo.m_spectator[spectators.size()],momentum,{},status));
	}

    MOCK vParticles Particles() const;
	vParticles allHadrons() const;
    crefParticles Protons(ParticleStatus = ParticleStatus::any) const;
    refParticles Protons(ParticleStatus = ParticleStatus::any);
    crefParticles Pions(ParticleStatus = ParticleStatus::any) const;
    refParticles Pions(ParticleStatus = ParticleStatus::any);
    crefParticles Neutrons(ParticleStatus = ParticleStatus::any) const;
    refParticles Neutrons(ParticleStatus = ParticleStatus::any);

    MOCK const EventHistory &History() const { return m_history; }
    EventHistory &History() { return m_history; }

    bool operator==(const Event &other) const {
        return m_nuc == other.m_nuc && m_remnant == other.m_remnant //&& m_mom == other.m_mom
				&&leptonsIn==other.leptonsIn && leptonsOut==other.leptonsOut
				&&hadronsIn==other.hadronsIn && hadronsOut==other.hadronsOut;
    }

    int &ProcessId() { return m_process_id; }
    const int &ProcessId() const { return m_process_id; }

  private:
    // Helper functions
    template <class UnaryPred>
    crefParticles FilterParticles(const vParticles &particles, UnaryPred pred) const {
        crefParticles result;
        std::copy_if(particles.begin(), particles.end(), std::back_inserter(result), pred);
        return result;
    }
    template <class UnaryPred> refParticles FilterParticles(vParticles &particles, UnaryPred pred) {
        refParticles result;
        std::copy_if(particles.begin(), particles.end(), std::back_inserter(result), pred);
        return result;
    }

    // Variables
	ProcessInfo& m_processInfo;
    std::shared_ptr<Nucleus> m_nuc;
    NuclearRemnant m_remnant{};
    //vMomentum m_mom{};
    double m_wgt{};
    vParticles leptonsIn{}, leptonsOut{}, hadronsIn{}, hadronsOut{}, spectators{};
    EventHistory m_history{};
    double flux{};
    int m_process_id{};
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