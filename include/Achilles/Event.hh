#ifndef EVENT_HH
#define EVENT_HH

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "Achilles/Achilles.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/NuclearRemnant.hh"
#include "Achilles/ProcessInfo.hh"

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
    Event(std::shared_ptr<Nucleus>, std::vector<FourVector>, double);
    Event(const Event &);
    Event &operator=(const Event &);
    MOCK ~Event() = default;

    void Finalize();

    MOCK const NuclearRemnant &Remnant() const { return m_remnant; }

    MOCK const vMomentum &Momentum() const { return m_mom; }
    MOCK vMomentum &Momentum() { return m_mom; }

    MOCK const std::shared_ptr<Nucleus> &CurrentNucleus() const { return m_nuc; }
    MOCK std::shared_ptr<Nucleus> &CurrentNucleus() { return m_nuc; }

    const double &Flux() const { return flux; }
    double &Flux() { return flux; }

    MOCK vParticles Particles() const;
    MOCK const vParticles &Hadrons() const { return m_hadrons; }
    MOCK vParticles &Hadrons() { return m_hadrons; }
    MOCK const vParticles &Leptons() const { return m_leptons; }
    MOCK vParticles &Leptons() { return m_leptons; }
    MOCK const double &Weight() const { return m_wgt; }
    MOCK double &Weight() { return m_wgt; }
    void Rotate(const std::array<double, 9> &);
    void Display() const;

    crefParticles Protons(ParticleStatus = ParticleStatus::any) const;
    refParticles Protons(ParticleStatus = ParticleStatus::any);
    crefParticles Neutrons(ParticleStatus = ParticleStatus::any) const;
    refParticles Neutrons(ParticleStatus = ParticleStatus::any);

    MOCK const EventHistory &History() const { return m_history; }
    EventHistory &History() { return m_history; }

    bool operator==(const Event &other) const {
        return m_nuc == other.m_nuc && m_remnant == other.m_remnant && m_mom == other.m_mom &&
               m_leptons == other.m_leptons;
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
    std::shared_ptr<Nucleus> m_nuc;
    NuclearRemnant m_remnant{};
    vMomentum m_mom{};
    double m_wgt{};
    vParticles m_leptons{}, m_hadrons{};
    EventHistory m_history{};
    double flux{};
    int m_process_id{};
};

} // namespace achilles

#endif
