#ifndef EVENT_HH
#define EVENT_HH

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "Achilles/Achilles.hh"
#include "Achilles/HardScatteringEnum.hh"
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

class Event {
  public:
    Event(double vWgt = 0) : m_wgt{vWgt} {}
    Event(std::shared_ptr<Nucleus>, std::vector<FourVector>, double);
    MOCK ~Event() = default;

    void SetHardScatteringType(HardScatteringType type) { m_type = type; }
    void Finalize();

    MOCK const NuclearRemnant &Remnant() const { return m_remnant; }

    MOCK const vMomentum &Momentum() const { return m_mom; }
    MOCK vMomentum &Momentum() { return m_mom; }

    const std::shared_ptr<Nucleus> &CurrentNucleus() const { return m_nuc; }
    MOCK std::shared_ptr<Nucleus> &CurrentNucleus() { return m_nuc; }

    const double &Flux() const { return flux; }
    double &Flux() { return flux; }

    MOCK vParticles Particles() const;
    const vParticles &Hadrons() const;
    MOCK vParticles &Hadrons();
    const vParticles &Leptons() const { return m_leptons; }
    vParticles &Leptons() { return m_leptons; }
    MOCK const double &Weight() const { return m_wgt; }
    MOCK double &Weight() { return m_wgt; }
    void Rotate(const std::array<double, 9> &);
    void Display() const;

    bool operator==(const Event &other) const {
        return m_type == other.m_type && m_nuc == other.m_nuc && m_remnant == other.m_remnant &&
               m_mom == other.m_mom && m_me == other.m_me && m_leptons == other.m_leptons;
    }

  private:
    HardScatteringType m_type{HardScatteringType::None};
    std::shared_ptr<Nucleus> m_nuc;
    NuclearRemnant m_remnant{};
    vMomentum m_mom{};
    std::vector<double> m_me;
    double m_wgt{};
    vParticles m_leptons{};
    vParticles m_history{};
    double flux;
};

} // namespace achilles

#endif
